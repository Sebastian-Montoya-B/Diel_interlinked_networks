
#------------------------------------------------------------------------------#
#                                                                              #
#                 1. Loading packages and sourcing functions                   #
#                                                                              #
#------------------------------------------------------------------------------#
# To install the infomapecology package for the first time follow https://infomapecology.ecomplab.com/installation.html

packs<-c("infomapecology",
         "emln",
         "igraph",
         "magrittr",
         "scales",
         "stringr",
         "tidyverse",
         "tnet", 
         "dendextend",
         "furrr",
         "purrr",
         "progressr",
         "betapart")

for (i in 1:length(packs)){
  
  if(!require(packs[i],character.only=TRUE)){
    install.packages(packs[i])
  }
  library(packs[i],character.only=TRUE)
}

source("Code/00_helping_functions.R")


df<-read.csv("Data/interactions_with_inttype.csv")
polli<-df

plist<-read.csv("Data/plant_abundances.csv", sep=",")

if (length(which(is.na(plist$abundance)))>0){
  plist<-plist[-which(is.na(plist$abundance)),]  
}

plist


df2<-plist[,
           c("ID", "date", "resource", "abundance")] %>%
  split(., .$ID)

df2 <-
  lapply(df2, function(x) {
    x %>%
      group_by(resource, date) %>%
      summarize(abundance = max(abundance, na.rm=TRUE), .groups = "drop")
  })

df2_clean <- lapply(df2, function(x) {
  x %>%
    group_by(resource) %>%
    summarize(abundance = sum(abundance, na.rm=TRUE), .groups = "drop")
})

#------------------------------------------------------------------------------#
#                                                                              #
#                               5. Plant check                                 #
#                                                                              #
#------------------------------------------------------------------------------#


polli<-polli[-which(polli$int_type=="non-trophic"),]
polli.div_raw<-polli %>% 
  split(f=polli$year) %>% 
  lapply(function(x){split(x, f=x$site_code)}) %>% 
  unlist(., recursive=FALSE) 



samp_intensity<-
  lapply(polli.div_raw, function(x){x$resource %>% 
      table() %>% 
      data.frame() %>% 
      `colnames<-`(c("resource", "intensity"))})%>%
  bind_rows(., .id="source") %>% group_by (source)%>%
  summarise( intensity= sum(intensity))

samp_intensity2<-
  lapply(polli.div_raw, function(x){x$resource %>% 
      table() %>% 
      data.frame() %>% 
      `colnames<-`(c("resource", "intensity"))})%>%
  bind_rows(., .id="source") %>% group_by (source) %>% 
  separate(., source, into=c("year", "site"), sep="\\.", extra="merge")

# plants that appear in the interactions but have no measured abundance
missing_plants <- map2(polli.div_raw, df2_clean,
                       ~ setdiff(unique(.x$resource), .y$resource))

imputed_share <- map2_dfr(polli.div_raw, missing_plants, ~ tibble(
  n_int            = nrow(.x),                    # one row = one interaction
  n_int_imputed    = sum(.x$resource %in% .y),
  n_plants         = n_distinct(.x$resource),
  n_plants_imputed = length(.y),
  plants_imputed   = paste(sort(.y), collapse = "; ")
), .id = "source") %>%
  mutate(prop_int    = n_int_imputed / n_int,
         prop_plants = n_plants_imputed / n_plants) %>%
  separate(source, into = c("year", "site"), sep = "\\.", extra = "merge")

imputed_share %>% arrange(desc(prop_int))

# pooled across all networks
imputed_share %>%
  summarise(prop_int_overall = sum(n_int_imputed) / sum(n_int),
            mean_prop_int    = mean(prop_int),
            max_prop_int     = max(prop_int))


df2_clean<- ## Give a value of minimum abundance for missing plants
  pmap(list(polli.div_raw, df2_clean), ~ 
         if (length(setdiff(..1$resource,..2$resource))>0){
           bind_rows(..2,
                     tibble(resource=setdiff(..1$resource,..2$resource),
                            abundance=min(..2$abundance)))
         } else{
           ..2
         }
  )




#------------------------------------------------------------------------------#
#                                                                              #
#                       5. Splitting data into layers                          #
#                                                                              #
#------------------------------------------------------------------------------#

polli.div_raw<-polli.div_raw %>%
  lapply(.,function(x) {split(x,f=x$record_period)}) %>%
  unlist(., recursive=FALSE) %>% 
  lapply(.,function(x) {split(x,f=x$int_type)})


polli.div_clean<-
  lapply(polli.div_raw,function(y){lapply(y,function(z){z[,c("consumer","resource")]%>% 
      group_by_all() %>% 
      count})}) %>%
  lapply(function(y){lapply(y,function(z){`colnames<-`(z,c("from", "to", "weight"))})}) 

str(polli.div_clean, max.level=1 )


samp_intensity3<- lapply(polli.div_clean, function (x) {
  br<-bind_rows(x, .id="int_type") 
  (sum(br$weight)) %>% 
    data.frame(intensity=.)
})%>% bind_rows(., .id="source") %>%
  separate(., source, into=c("year", "site", "period"), sep="\\.",extra="merge")


per_ant<- lapply(polli.div_clean, function (x) {
  br<-bind_rows(x, .id="int_type") %>%
    group_by(., int_type)%>%
    summarize(., intensity=sum(weight))%>%
    data.frame(.)
}) %>% bind_rows(., .id="source")%>%
  separate(., source, into=c("year", "site", "period"), sep="\\.",extra="merge")



polli.div<-polli.div_clean %>%
  lapply(.,function(y){lapply(y,
                              function(z){z$weight<-z$weight/sum(z$weight)
                              return(z) })})



#------------------------------------------------------------------------------#
#                                                                              #
#                       5. Creating layer attributes                           #
#                                                                              #
#------------------------------------------------------------------------------#

layer_attrib <-
  polli.div %>% 
  lapply(., 
         function(y){
           tibble(layer_id=1:length(names(y)),layer_name=names(y))
         }
  )


#------------------------------------------------------------------------------#
#                                                                              #
#                       6. Creating interlayer links                           #
#                                                                              #
#------------------------------------------------------------------------------#

#### INTERLAYER LINKS FOR PLANTS

df1<-polli.div %>%
  lapply(., function(x) do.call(rbind,x) %>%
           .$to %>% 
           unique %>% 
           data.frame(resource=.))

# Extract prefixes from df1 names (correctly capturing up to the second period)
prefixes_df1 <- gsub("^([^.]+\\.[^.]+).*", "\\1", names(df1))

# Split df1 into a list based on the prefixes
df1_split <- split(df1, prefixes_df1) %>% lapply(.,length) %>% unlist

df2_clean_rep<- rep(df2_clean, times=df1_split)



fabun<-
  pmap(list(df1,df2_clean_rep),~..1 %>%
         left_join(..2, by = "resource")
)

fabun2<-
  fabun %>% lapply(., function(x) {x$abundance<-x$abundance/sum(x$abundance)
  return(x)}) 


interlayer<-
  pmap(list(fabun2, layer_attrib), ~ 
         data.frame(layer_from=rep(..2$layer_name, each=nrow(..1)),
                    node_from=rep(..1$resource, times=nrow(..2)),
                    layer_to=rep(..2$layer_name %>% reorder_vector(.), each=nrow(..1)),
                    node_to=rep(..1$resource, times=nrow(..2)),
                    weight=rep(..1$abundance, times=nrow(..2))))




#------------------------------------------------------------------------------#
#                                                                              #
#                      7. Creating Multilayer networks                         #
#                                                                              #
#------------------------------------------------------------------------------#

multi.polli<-pmap(list(polli.div, interlayer,layer_attrib), 
                  ~create_multilayer_network(list_of_layers = ..1,
                                             layer_attributes = ..3, 
                                             interlayer_links = ..2,
                                             bipartite = TRUE,
                                             directed = FALSE
                  ))


#------------------------------------------------------------------------------#
#                                                                              #
#                      8. Exporting processed objects                          #
#                                                                              #
#------------------------------------------------------------------------------#

# Save specie slist
#saveRDS(df$resource %>% unique, "Data/Processed objects/plant_list.rds") # Original code using the whole raw data
saveRDS(df$resource %>% unique, "Data/Processed objects/subset/plant_list.rds")


# Save visitor abundances
insect_ab<-
  polli.div_clean %>% 
  lapply(., function(x){
    lapply(x, function(y){ y %>%
        group_by(from) %>%
        summarize(insect_abundance=sum(weight))
      
    })})


# Save sampling intensity data per site and period
#saveRDS(samp_intensity3, "Data/Processed objects/samp_intensity3.rds") # Original code using the whole raw data
saveRDS(samp_intensity3, "Data/Processed objects/subset/samp_intensity3.rds")


igraph_layers<-
  lapply(multi.polli, function(x) get_igraph(x, bipartite=TRUE, directed=FALSE)$layers_igraph)


#------------------------------------------------------------------------------#
#                                                                              #
#                          9. Preparing data for Figure 2                      #
#                                                                              #
#------------------------------------------------------------------------------#

iwint<-df

## How many visitors do we have from each order
### mutualists
table(iwint$con_order[
  iwint$int_type=="mutualistic"] )
### antagonists
table(iwint$con_order[
  iwint$int_type=="antagonistic"] )


iwint %>% dplyr::select(year, site_code, record_period, int_type, consumer)->smdf
smdf<-smdf%>% unique

smdf %>% group_by(year, site_code, record_period) %>%
  summarise(mutualistic = sum(int_type=="mutualistic"),
            antagonistic = sum(int_type=="antagonistic"))->smdf

smdf$source<-paste0(smdf$year, ".", smdf$site_code)

#saveRDS(smdf,"Data/Processed objects/Fig_2.rds") # Original code using the whole raw data
saveRDS(smdf,"Data/Processed objects/subset/Fig_2.rds")




#------------------------------------------------------------------------------#
#                                                                              #
#                         10. Preparing data for Figure 3                      #
#                                                                              #
#------------------------------------------------------------------------------#


####### Estimate unweighted degree

layer_degree <- igraph_layers %>% 
  lapply(., function(x) {
    lapply(names(x), function(y_name) { 
      y <- x[[y_name]]  # Access the actual element in the list
      
      # Convert y to a biadjacency matrix
      mat <- as_biadjacency_matrix(y)
      
      # Ensure mat is numeric before performing any calculations
      if (!is.numeric(mat)) {
        mat <- as.numeric(mat)  # Convert to numeric if necessary
      }
      
      # Check if the matrix is empty or malformed (i.e., has no columns or rows)
      if (ncol(mat) == 0 || nrow(mat) == 0) {
        stop("The biadjacency matrix has no rows or columns!")
      }
      
      # Normalize by the number of columns (if ncol > 0)
      normalized_result <- rowSums(mat)
      
      # Convert the result to a data frame, set column name and row names
      result_df <- as.data.frame(normalized_result) %>%
        `colnames<-`(y_name) %>%
        rownames_to_column(var = "resource")
      
      return(result_df)
    })
  })




degree_summ<-list()
for (i in 1:length(layer_degree)){
  degree_summ[[i]]<-
    Reduce(function(x,y) {
      merge(x, y, by.x="resource", by.y="resource",
            all=TRUE)
    },layer_degree[[i]])
  degree_summ[[i]][is.na(degree_summ[[i]])]<-0
}
names(degree_summ)<-names(layer_degree)
degree_summ


degree_summ<-pmap(list(degree_summ, df2_clean_rep), ~merge(..1,..2, by="resource"))



ant_vs_mut_cor <- calculate_partial_cor(degree_summ)

size<-
  sapply(degree_summ, nrow) %>% data.frame(source=names(.), size=.) %>%
  separate(., source, into=c("year", "site", "period"), sep="\\.", extra="merge")

#

# Combine the list of data frames into a single data frame
ant_vs_mut_cor_df <- ant_vs_mut_cor %>%
  separate(name, into = c("year", "site", "period"), sep = "\\.", extra = "merge")%>%
  merge(., size, by=c("year", "site", "period")) %>%
  merge(., samp_intensity3, by=c("year", "site", "period"))



sil<-
  lapply(seq_along(degree_summ), function(i) {
    
    # Extract the relevant columns
    data <- cbind(
      antagonistic = degree_summ[[i]]$antagonistic,
      mutualistic = degree_summ[[i]]$mutualistic,
      abundance = (degree_summ[[i]]$abundance)
    ) %>% as.data.frame(.)
  })
names(sil)<-names(degree_summ)

set.seed(14)
null_sil <- lapply(sil, function(x) {
  tryCatch({
    
    generate_null_models(x, 
                         columns_to_shuffle = c("antagonistic", "mutualistic"), 
                         num_models = 1000)
  }, 
  warning = function(w) {
    message("Warning for element: ", w)
    return(NULL)  # or handle warnings differently
  },
  error = function(e) {
    message("Error for element: ", e$message)
    return(NULL)  # Skip problematic elements
  })
})



ant_vs_mut_cor_null <- lapply(null_sil, function(network) {
  lapply(network, function(df) {
    calculate_partial_cor(
      list(temp = df),  # Wrap single DF in a list
      ant_col = "antagonistic",
      mut_col = "mutualistic",
      abund_col = "abundance",
      verbose = FALSE
    )
  })%>% bind_rows(.)
}) 


ant_vs_mut_cor_null$`2014.A.1.morning`
ant_vs_mut_cor_df


z_score<-list()
null_mean<-list()
effect_size<-list()
for(i in 1:length(ant_vs_mut_cor_null)){
  if (nrow(ant_vs_mut_cor_null[[i]])>0 & !is.na(ant_vs_mut_cor_df$correlation[i])){
    
    
    # Calculate the mean of the null distribution
    null_mean[[i]] <- mean(na.omit(ant_vs_mut_cor_null[[i]]$correlation))
    # Z-score
    z_score[[i]] <- ((ant_vs_mut_cor_df$correlation[i]) - mean(na.omit(ant_vs_mut_cor_null[[i]]$correlation))) / sd(na.omit(ant_vs_mut_cor_null[[i]]$correlation))
    # Standardized effect size
    effect_size[[i]] <- (ant_vs_mut_cor_df$correlation[i]) - mean(na.omit(ant_vs_mut_cor_null[[i]]$correlation))
  } else {
    z_score[[i]]<-NA
    effect_size[[i]]<-NA
    null_mean[[i]]<-NA
  }
}

ant_vs_mut_cor_df$z_score<-z_score %>% do.call(c,.)
ant_vs_mut_cor_df$null_mean<-null_mean %>% do.call(c,.)
ant_vs_mut_cor_df$effect_size <-effect_size %>% do.call(c,.)
ant_vs_mut_cor_df<-ant_vs_mut_cor_df[which(ant_vs_mut_cor_df$status=="success"),] # Select only those networks with successful correlation estimation


ant_vs_mut_cor_df %>%  dplyr::select(year, site, period, correlation, n_species,
                                     size,  z_score,null_mean, effect_size)



full_data<- ant_vs_mut_cor_df
full_data$ID<-seq(1, nrow(full_data))

#saveRDS(full_data, "Data/Processed objects/degree_unweighted.rds") # Original code using the whole raw data
saveRDS(full_data, "Data/Processed objects/subset/degree_unweighted.rds") 


####### Estimate weighted degree

layer_as_matrix2 <- igraph_layers %>%
  lapply(., function(x) {
    lapply(names(x), function(y_name) { 
      y <- x[[y_name]]  # Access the actual element in the list
      
      my<-as_weighted_biadjacency_matrix(y) %>%t()
      names(my)<-names(y)
      return(my[,order(colnames(my)), drop=FALSE])
    })
  }) 




for (i in 1:length(layer_as_matrix2)){
  names(layer_as_matrix2[[i]])<-names(igraph_layers[[i]])
}

layer_as_matrix2

layer_pref_deg <- pmap(list(layer_as_matrix2, insect_ab), ~ 
                         pmap(list(..1, ..2), ~ {
                           pref_Shannon(..1, ..2$insect_abundance) %>%
                             as.data.frame() %>%
                             rownames_to_column(var = "resource")
                         })
)



pref_deg_summ<-list()
for (i in 1:length(layer_pref_deg)){
  pref_deg_summ[[i]]<-
    Reduce(function(x,y) {
      merge(x, y, by.x="resource", by.y="resource",
            all=TRUE)
    },layer_pref_deg[[i]])
  pref_deg_summ[[i]][is.na(pref_deg_summ[[i]])]<-0
}


names(pref_deg_summ)<-names(layer_as_matrix2)

for (i in 1:length(pref_deg_summ)){
  colnames(pref_deg_summ[[i]])[seq(2, ncol(pref_deg_summ[[i]]))]<-names(igraph_layers[[i]])
}



pref_deg_summ<-pmap(list(pref_deg_summ, df2_clean_rep), ~merge(..1,..2, by="resource"))

if (any(which(sapply(pref_deg_summ, function(x) (nrow(x)))<3))){
  pref_deg_summ<-pref_deg_summ[-which(sapply(pref_deg_summ, function(x) (nrow(x)))<3)]  
}



ant_vs_mut_cor2 <- calculate_partial_cor(pref_deg_summ)


ant_vs_mut_cor_df2 <- ant_vs_mut_cor2%>%
  separate(name, into = c("year", "site", "period"), sep = "\\.", extra = "merge")%>%
  merge(., size, by=c("year", "site", "period")) %>%
  merge(., samp_intensity3, by=c("year", "site", "period"))




sil2<-
  lapply(seq_along(pref_deg_summ), function(i) {
    
    # Extract the relevant columns
    data <- cbind(
      antagonistic = pref_deg_summ[[i]]$antagonistic,
      mutualistic = pref_deg_summ[[i]]$mutualistic,
      abundance = (pref_deg_summ[[i]]$abundance)
    ) %>% as.data.frame(.)
  })
names(sil2)<-names(pref_deg_summ)

set.seed(14)
null_sil2 <- lapply(sil2, function(x) {
  tryCatch({
    
    generate_null_models(x, 
                         columns_to_shuffle = c("antagonistic", "mutualistic"), 
                         num_models = 1000)
  }, 
  warning = function(w) {
    message("Warning for element: ", w)
    return(NULL)  # or handle warnings differently
  },
  error = function(e) {
    message("Error for element: ", e$message)
    return(NULL)  # Skip problematic elements
  })
})


ant_vs_mut_cor_null2 <- lapply(null_sil2, function(network) {
  lapply(network, function(df) {
    calculate_partial_cor(
      list(temp = df),  # Wrap single DF in a list
      ant_col = "antagonistic",
      mut_col = "mutualistic",
      abund_col = "abundance",
      verbose = FALSE
    ) 
  })%>% bind_rows(.)
}) 


ant_vs_mut_cor_null2$`2014.A.1.morning`
ant_vs_mut_cor_df2

z_score2<-list()
null_mean2<-list()
effect_size2<-list()
for(i in 1:length(ant_vs_mut_cor_null2)){
  if (nrow(ant_vs_mut_cor_null2[[i]])>0 & !is.na(ant_vs_mut_cor_df2$correlation[i])){
    
    # Calculate the mean of the null distribution
    null_mean2[[i]] <- mean(na.omit(ant_vs_mut_cor_null2[[i]]$correlation))
    # Z-score
    z_score2[[i]] <- ((ant_vs_mut_cor_df2$correlation[i]) - mean(na.omit(ant_vs_mut_cor_null2[[i]]$correlation))) / sd(na.omit(ant_vs_mut_cor_null2[[i]]$correlation))
    # Standardized effect size
    effect_size2[[i]] <- (ant_vs_mut_cor_df2$correlation[i]) - mean(na.omit(ant_vs_mut_cor_null2[[i]]$correlation))
  } else {
    z_score2[[i]]<-NA
    null_mean2[[i]]<-NA
    effect_size2[[i]]<-NA
  }
}

ant_vs_mut_cor_df2$z_score<-z_score2 %>% do.call(c,.)
ant_vs_mut_cor_df2$null_mean<-null_mean2 %>% do.call(c,.)
ant_vs_mut_cor_df2$effect_size <-effect_size2 %>% do.call(c,.)
ant_vs_mut_cor_df2<-ant_vs_mut_cor_df2[which(ant_vs_mut_cor_df2$status=="success"),]



ant_vs_mut_cor_df2



full_data2<- ant_vs_mut_cor_df2
full_data2$ID<-seq(1, nrow(full_data2))



full_data2$source<-paste0(full_data2$year, full_data2$site)

#saveRDS(full_data2, "Data/Processed objects/degree_weighted.rds") # Original code using the whole raw data
saveRDS(full_data2, "Data/Processed objects/subset/degree_weighted.rds") 

#------------------------------------------------------------------------------#
#                                                                              #
#              10. Generating null models for modularity analysis              #
#                                                                              #
#------------------------------------------------------------------------------#

# This part may take some minutes

AM_mat<-polli.div_clean%>%
  map_depth(2, edge_list_to_bipartite)

str(AM_mat, max.level=1)

AM_mat$`2014.A.1.morning`

# Creating a function to generate null models based on Patefield's algorithm (r2d method in bipartite package)
create_null_models <- function(N = 1000) {
  AM_mat_null <- AM_mat %>%
    map_depth(2, ~ {
      tryCatch({
        
        null_matrices <- bipartite::nullmodel(.x, method = "r2d", N = N)
        
        # Copy row and column names from original to all null matrices
        null_matrices <- lapply(null_matrices, function(mat) {
          dimnames(mat) <- dimnames(.x)
          return(mat)
        })
        
        return(null_matrices)
        
      }, error = function(e) {
        # Return original matrix replicated N times (with names already intact)
        rep(list(.x), N)
      })
    })
  return(AM_mat_null)
}

# Usage
set.seed(265)
AM_mat_null <- create_null_models(N = 1000)
null_layers<-AM_mat_null %>% map_depth(3, ~ {return(.x/sum(.x))})

# Checking a resulting null network
null_layers$`2014.A.1.morning`$antagonistic[[1]]
null_layers$`2014.A.1.morning`$mutualistic[[1]]



# Building multilayer networks for null models in parallel
handlers(global = TRUE)
handlers("progress")

plan(multisession, workers = parallel::detectCores() - 1)

with_progress({
  p <- progressor(along = null_layers)
  
  null_multi <- future_imap(null_layers, function(layer_set, u) {
    p()  # Update progress
    future_map(1:length(layer_set$antagonistic), function(i) {
      lay_pairs <- list(
        antagonistic = layer_set$antagonistic[[i]], 
        mutualistic = layer_set$mutualistic[[i]]
      )
      
      create_multilayer_network(
        list_of_layers = lay_pairs,
        layer_attributes = layer_attrib[[u]], 
        interlayer_links = interlayer[[u]],
        bipartite = TRUE,
        directed = FALSE
      )
    }, .options = furrr_options(seed = TRUE))
  }, .options = furrr_options(seed = TRUE))
})


safe_run_infomap <- function(network) {
  tryCatch({
    run_infomap_multilayer(network, flow_model = "undirected")
  }, error = function(e) {
    message("Infomap failed: ", e$message)
    return(NULL)
  })
}

# Add progress with purrr's progress option
set.seed(15)
map_null <- null_multi %>% 
  map_depth(2, safe_run_infomap, .progress = TRUE)


#saveRDS(map_null,"Data/Processed objects/map_null.rds") # Original code using the whole raw data
saveRDS(map_null,"Data/Processed objects/subset/map_null.rds") 


#------------------------------------------------------------------------------#
#                                                                              #
#                         11. Calculate observed modularity                    #
#                                                                              #
#------------------------------------------------------------------------------#

set.seed(24)
map_all<-lapply(multi.polli,run_infomap_multilayer, flow_model="undirected")
#saveRDS(map_all,"Data/Processed objects/map_all.rds") # Original code using the whole raw data
saveRDS(map_all,"Data/Processed objects/subset/map_all.rds") 


#------------------------------------------------------------------------------#
#                                                                              #
#                    12. Calculate beta diversity of plants                    #
#                                                                              #
#------------------------------------------------------------------------------#



mas<-
  igraph_layers %>%
  map_depth(., 2, ~ .x %>% 
              as_biadjacency_matrix,attr = "weight")

combine_biadjacency <- function(mat_list) {
  mat_list |>
    map(~ as_tibble(.x, rownames = "row_id")) |>   # matrix -> tibble, keep rownames
    reduce(full_join, by = "row_id") |>            # union of row names, NA where absent
    mutate(across(-row_id, ~ replace_na(.x, 0))) |> # fill missing entries with 0
    column_to_rownames("row_id") |>
    as.matrix()
}

aggregate_mas<-lapply(mas,function(x) combine_biadjacency(x))
str(aggregate_mas)

nested_mas <- aggregate_mas |>
  split(str_extract(names(aggregate_mas), "^[0-9]{4}"))
nested_mas <- nested_mas |>
  map(\(yr_list) split(yr_list, str_split_i(names(yr_list), "\\.", 2)))
nested_mas %>% str

beta_mas<-nested_mas[c(1,3)]


plant_beta <- function(webs) {
  # plants interacting in each period (union of both layers, from combine_biadjacency)
  plants <- map(webs, \(m) rownames(m)[rowSums(m) > 0])
  all_p  <- unique(unlist(plants))
  
  pa <- t(sapply(plants, \(p) as.integer(all_p %in% p)))
  dimnames(pa) <- list(names(webs), all_p)
  
  bp    <- beta.pair(pa, index.family = "sorensen")
  pairs <- combn(names(webs), 2)          # same order as as.vector(dist)
  
  tibble(
    i        = pairs[1, ],
    j        = pairs[2, ],
    n_i      = lengths(plants)[pairs[1, ]],
    n_j      = lengths(plants)[pairs[2, ]],
    n_shared = map2_int(pairs[1, ], pairs[2, ],
                        \(a, b) length(intersect(plants[[a]], plants[[b]]))),
    sor      = as.vector(bp$beta.sor),    # total plant dissimilarity
    sim      = as.vector(bp$beta.sim),    # turnover (replacement)
    sne      = as.vector(bp$beta.sne)     # nestedness (loss of a subset)
  )
}

plant_beta_out <- map_depth(beta_mas, 2, plant_beta)

plant_beta_df <- plant_beta_out |>
  map(bind_rows, .id = "site") |>
  bind_rows(.id = "year") |>
  mutate(
    comparison = case_when(
      str_detect(i, "morning")   & str_detect(j, "afternoon") ~ "morning_afternoon",
      str_detect(i, "morning")   & str_detect(j, "night")     ~ "morning_night",
      str_detect(i, "afternoon") & str_detect(j, "night")     ~ "afternoon_night",
      TRUE ~ NA_character_
    ),
    sne_prop = sne / sor                  # share of dissimilarity due to nestedness
  )

print(plant_beta_df, n = Inf)

# summary per comparison type
plant_beta_df |>
  group_by(comparison) |>
  summarise(across(c(n_i, n_j, n_shared, sor, sim, sne, sne_prop),
                   \(x) mean(x, na.rm = TRUE)))

#------------------------------------------------------------------------------#
#                                                                              #
#   13. Which plants keep interacting at night, and do they keep their role?   #
#                                                                              #
#------------------------------------------------------------------------------#

# 
#
# Analysis 1 - WHO is shared?
#   Are the plants that also interact at night ("shared") the afternoon plants
#   with most consumers (and/or the most abundant plants), or a random subset?
#   Per site: AUC = probability that a randomly chosen shared plant has a higher
#   afternoon value than a randomly chosen afternoon-only plant (ties = 0.5).
#     AUC = 0.5 -> shared plants are a random subset of afternoon plants
#     AUC > 0.5 -> night keeps the afternoon "hubs" (well-visited / abundant)
#     AUC < 0.5 -> night keeps the poorly visited / rare plants
#   Across sites: mean AUC tested with a permutation test that shuffles the
#   shared / afternoon-only labels WITHIN each site (keeps each site's number
#   of shared plants fixed).
#   Also reported: share of afternoon interactions held by shared plants vs
#   their share of afternoon species.
#
# Analysis 2 - do shared plants keep their ROLE?
#   Among shared plants only, Spearman correlation between afternoon and night
#   degree, separately for the mutualistic and the antagonistic layer.
#     High rho -> plants that are hubs in the afternoon are also hubs at night
#     Low rho  -> plants re-rank: consumers use the same plants differently
#   Sites with < 4 shared plants are skipped (rho not meaningful).
#

# =============================================================================

site_of   <- function(nm) sub("\\.[0-9]+\\.(morning|afternoon|night)$", "", nm)
period_of <- function(nm) sub("^.*\\.[0-9]+\\.", "", nm)
active    <- function(df) df[(df$antagonistic + df$mutualistic) > 0, ]

# AUC (Mann-Whitney U / (n1 * n2)); ties count 0.5
auc <- function(x_shared, x_other) {
  if (length(x_shared) == 0 || length(x_other) == 0) return(NA_real_)
  mean(outer(x_shared, x_other, function(a, b) (a > b) + 0.5 * (a == b)))
}

# Afternoon / night pair for each site that has both
site_pairs <- function(degree_summ) {
  nm <- names(degree_summ)
  out <- list()
  for (s in unique(site_of(nm))) {
    aft <- nm[site_of(nm) == s & period_of(nm) == "afternoon"]
    nig <- nm[site_of(nm) == s & period_of(nm) == "night"]
    if (length(aft) == 1 && length(nig) == 1)
      out[[s]] <- list(aft = active(degree_summ[[aft]]),
                       nig = active(degree_summ[[nig]]))
  }
  out
}

# ---- Analysis 1: who is shared? ---------------------------------------------
shared_identity <- function(degree_summ) {
  pairs <- site_pairs(degree_summ)
  rows <- lapply(names(pairs), function(s) {
    A <- pairs[[s]]$aft; N <- pairs[[s]]$nig
    A$shared <- A$resource %in% N$resource
    A$total  <- A$antagonistic + A$mutualistic
    sh <- A[A$shared, ]; ao <- A[!A$shared, ]
    data.frame(
      site = s,
      n_afternoon = nrow(A), n_shared = nrow(sh), n_afternoon_only = nrow(ao),
      share_species      = nrow(sh) / nrow(A),
      share_int_total    = sum(sh$total)        / sum(A$total),
      share_int_mut      = sum(sh$mutualistic)  / max(sum(A$mutualistic), 1),
      share_int_ant      = sum(sh$antagonistic) / max(sum(A$antagonistic), 1),
      auc_total     = auc(sh$total,        ao$total),
      auc_mut       = auc(sh$mutualistic,  ao$mutualistic),
      auc_ant       = auc(sh$antagonistic, ao$antagonistic),
      auc_abundance = auc(sh$abundance,    ao$abundance))
  })
  do.call(rbind, rows)
}

# Stratified permutation test: is the mean AUC across sites different from
# what random within-site subsets of the same size would give?
perm_test_auc <- function(degree_summ, variable = c("total", "mutualistic",
                                                    "antagonistic", "abundance"), n_perm = 9999, seed = 1) {
  variable <- match.arg(variable)
  pairs <- site_pairs(degree_summ)
  dat <- lapply(pairs, function(p) {
    A <- p$aft
    x <- if (variable == "total") A$antagonistic + A$mutualistic else A[[variable]]
    list(x = x, shared = A$resource %in% p$nig$resource)
  })
  dat <- dat[sapply(dat, function(d) any(d$shared) && any(!d$shared))]
  stat <- function(dl) mean(sapply(dl, function(d) auc(d$x[d$shared], d$x[!d$shared])))
  obs <- stat(dat)
  set.seed(seed)
  null <- replicate(n_perm, stat(lapply(dat, function(d) {
    d$shared <- sample(d$shared); d })))
  p <- (sum(abs(null - 0.5) >= abs(obs - 0.5)) + 1) / (n_perm + 1)   # two-sided
  data.frame(variable = variable, n_sites = length(dat),
             mean_auc = obs, null_mean = mean(null),
             null_2.5 = quantile(null, 0.025), null_97.5 = quantile(null, 0.975),
             p_value = p, row.names = NULL)
}

# ---- Analysis 2: do shared plants keep their role? --------------------------
role_persistence <- function(degree_summ, min_shared = 4) {
  pairs <- site_pairs(degree_summ)
  rows <- lapply(names(pairs), function(s) {
    A <- pairs[[s]]$aft; N <- pairs[[s]]$nig
    sh <- intersect(A$resource, N$resource)
    A <- A[match(sh, A$resource), ]; N <- N[match(sh, N$resource), ]
    r <- function(x, y) {
      if (length(sh) < min_shared || sd(x) == 0 || sd(y) == 0) return(NA_real_)
      suppressWarnings(cor(x, y, method = "spearman"))
    }
    data.frame(site = s, n_shared = length(sh),
               rho_mut = r(A$mutualistic,  N$mutualistic),
               rho_ant = r(A$antagonistic, N$antagonistic))
  })
  do.call(rbind, rows)
}

summarise_role <- function(rp) {
  f <- function(v) {
    x <- rp[[v]][!is.na(rp[[v]])]
    w <- suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE))
    data.frame(layer = sub("rho_", "", v), n_sites = length(x),
               mean_rho = mean(x), sd_rho = sd(x),
               n_positive = sum(x > 0), p_wilcoxon_vs_0 = w$p.value)
  }
  rbind(f("rho_mut"), f("rho_ant"))
}

# =============================================================================
# USAGE
# -----------------------------------------------------------------------------
id <- shared_identity(degree_summ); id
colMeans(id[, -1])                                   # quick overview
rbind(perm_test_auc(degree_summ, "total"),
      perm_test_auc(degree_summ, "mutualistic"),
      perm_test_auc(degree_summ, "antagonistic"),
      perm_test_auc(degree_summ, "abundance"))

rp <- role_persistence(degree_summ); rp
summarise_role(rp)

#Is the difference between the pollinator and the herbivore ranking significant?
wilcox.test(rp$rho_mut, rp$rho_ant, paired = TRUE, exact = FALSE)
# =============================================================================
#------------------------------------------------------------------------------#
#                                                                              #
#     14. Changes in the ranking of plants between afternoon and night         #
#                                                                              #
#------------------------------------------------------------------------------#

site_of   <- function(nm) sub("\\.[0-9]+\\.(morning|afternoon|night)$", "", nm)
period_of <- function(nm) sub("^.*\\.[0-9]+\\.", "", nm)

# ---- rank mismatch of every plant in one network ----------------------------
partial_mismatch <- function(ant, mut, abund, net = "") {
  n <- length(ant)
  if (n < 4) return(rep(NA_real_, n))          # too few plants for a partial correlation
  ra <- rank(ant); rm <- rank(mut); rb <- rank(abund)
  ea <- resid(lm(ra ~ rb)); em <- resid(lm(rm ~ rb))
  if (sd(ea) == 0 || sd(em) == 0) return(rep(NA_real_, n))
  d <- as.numeric(scale(ea) - scale(em))
  
  # check against ppcor: partial Spearman = 1 - sum(d^2) / (2(n-1))
  r_d     <- 1 - sum(d^2) / (2 * (n - 1))
  d
}




# ---- afternoon vs night for every shared plant ------------------------------
rank_shift <- function(degree_summ) {
  nm <- names(degree_summ)
  per_net <- function(df, net) {
    df <- df[(df$antagonistic + df$mutualistic) > 0, ]   # plants in the network
    data.frame(plant = df$resource,
               d = partial_mismatch(df$antagonistic, df$mutualistic,
                                    df$abundance, net))
  }
  rows <- list()
  for (s in unique(site_of(nm))) {
    aft <- nm[site_of(nm) == s & period_of(nm) == "afternoon"]
    nig <- nm[site_of(nm) == s & period_of(nm) == "night"]
    if (length(aft) != 1 || length(nig) != 1) next
    m <- merge(per_net(degree_summ[[aft]], aft), per_net(degree_summ[[nig]], nig),
               by = "plant", suffixes = c("_aft", "_night"))     # shared plants only
    if (nrow(m)) rows[[s]] <- cbind(community = s, m)
  }
  long <- do.call(rbind, rows)
  long$shift <- long$d_night - long$d_aft
  
  out <- do.call(rbind, lapply(split(long, long$plant), function(x) data.frame(
    plant             = x$plant[1],
    n_comm_shared     = nrow(x),                    # communities where shared
    n_comm_valid      = sum(!is.na(x$shift)),       # of those, with a valid d
    mean_d_aft        = mean(x$d_aft,   na.rm = TRUE),
    mean_d_night      = mean(x$d_night, na.rm = TRUE),
    mean_shift        = mean(x$shift,   na.rm = TRUE),
    n_herb_side_night = sum(x$d_night > 0, na.rm = TRUE),
    communities       = paste(x$community, collapse = ", "))))
  out$mean_shift[is.nan(out$mean_shift)] <- NA
  out <- out[order(-out$mean_shift), ]; rownames(out) <- NULL
  attr(out, "long") <- long
  out
}

# ---- visitors of selected plants --------------------------------------------
plant_visitors <- function(igraph_layers, plants) {
  nm <- names(igraph_layers)
  keep <- nm[period_of(nm) %in% c("afternoon", "night")]
  rows <- list()
  for (net in keep) for (layer in names(igraph_layers[[net]])) {
    g <- igraph_layers[[net]][[layer]]
    for (p in intersect(plants, V(g)$name)) {
      nb <- neighbors(g, p)$name
      if (length(nb))
        rows[[length(rows) + 1]] <- data.frame(
          plant = p, community = site_of(net), period = period_of(net),
          layer = layer, visitor = nb)
    }
  }
  do.call(rbind, rows)
}

rs <- rank_shift(degree_summ) 

#saveRDS(rs, "Data/Processed objects/rank_shift.rds") # Original code using the whole raw data
saveRDS(rs, "Data/Processed objects/subset/rank_shift.rds")


