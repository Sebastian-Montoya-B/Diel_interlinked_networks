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
         "progressr")

for (i in 1:length(packs)){
  
  if(!require(packs[i],character.only=TRUE)){
    install.packages(packs[i])
  }
  library(packs[i],character.only=TRUE)
}

source("Code/00_helping_functions.R")


df<-read.csv("Data/interactions_dualint.csv")
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


polli<-polli[-which(polli$Int_type2=="non-trophic"),]
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
  lapply(.,function(x) {split(x,f=x$Int_type2)})


polli.div_clean<-
  lapply(polli.div_raw,function(y){lapply(y,function(z){z[,c("consumer","resource")]%>% 
      group_by_all() %>% 
      count})}) %>%
  lapply(function(y){lapply(y,function(z){`colnames<-`(z,c("from", "to", "weight"))})}) 

str(polli.div_clean, max.level=1 )

#dividing dual interactions into antagonistic and mutualistic layers
empty_edges <- tibble(from = character(), to = character(), weight = numeric())

polli.div_merged <- map(polli.div_clean, function(x) { 
  
  if (is.null(x$dual)) return(x)
  
  dual_half <- x$dual %>%
    ungroup() %>%
    mutate(weight = weight / 2)
  
  merge_layer <- function(layer) {
    layer <- layer %||% empty_edges       
    layer %>%
      ungroup() %>%
      bind_rows(dual_half) %>%
      group_by(from, to) %>%
      summarise(weight = sum(weight), .groups = "drop")
  }
  
  list(
    antagonistic = merge_layer(x$antagonistic),
    mutualistic  = merge_layer(x$mutualistic)
  )
})




samp_intensity3<- lapply(polli.div_merged, function (x) {
  br<-bind_rows(x, .id="int_type") 
  (sum(br$weight)) %>% 
    data.frame(intensity=.)
})%>% bind_rows(., .id="source") %>%
  separate(., source, into=c("year", "site", "period"), sep="\\.",extra="merge")


per_ant<- lapply(polli.div_merged, function (x) {
  br<-bind_rows(x, .id="int_type") %>%
    group_by(., int_type)%>%
    summarize(., intensity=sum(weight))%>%
    data.frame(.)
}) %>% bind_rows(., .id="source")%>%
  separate(., source, into=c("year", "site", "period"), sep="\\.",extra="merge")



polli.div<-polli.div_merged %>%
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


make_interlayer_dual <- function(period, halve_dual = FALSE) {
  
  if (is.null(period$dual) || nrow(period$dual) == 0) return(NULL)
  
  empty_edges <- tibble::tibble(from = character(), weight = numeric())
  
  all_edges <- dplyr::bind_rows(
    dplyr::ungroup(period$antagonistic %||% empty_edges),
    dplyr::ungroup(period$mutualistic  %||% empty_edges),
    dplyr::ungroup(period$dual) %>%
      dplyr::mutate(weight = if (halve_dual) weight / 2 else weight)
  ) %>%
    dplyr::select(from, weight)
  
  community_total <- sum(all_edges$weight)
  
  rel <- all_edges %>%
    dplyr::filter(from %in% unique(period$dual$from)) %>%
    dplyr::group_by(from) %>%
    dplyr::summarise(w = sum(weight), .groups = "drop") %>%
    dplyr::mutate(weight = w / community_total) %>%
    dplyr::select(node = from, weight)
  
  dplyr::bind_rows(
    tibble::tibble(layer_from = "antagonistic", node_from = rel$node,
                   layer_to   = "mutualistic",  node_to   = rel$node,
                   weight = rel$weight),
    tibble::tibble(layer_from = "mutualistic",  node_from = rel$node,
                   layer_to   = "antagonistic", node_to   = rel$node,
                   weight = rel$weight)
  ) %>%
    as.data.frame()
}

insect.interlayer <- polli.div_clean %>%
  map(make_interlayer_dual) %>%
  compact()
###interlayer full

all_periods <- union(names(interlayer), names(insect.interlayer))

interlayer_full <- set_names(all_periods) %>%
  map(function(p) {
    bind_rows(interlayer[[p]], insect.interlayer[[p]])
  })

#------------------------------------------------------------------------------#
#                                                                              #
#                      7. Creating Multilayer networks                         #
#                                                                              #
#------------------------------------------------------------------------------#

multi.polli<-pmap(list(polli.div, interlayer_full,layer_attrib), 
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
#saveRDS(df$resource %>% unique, "Data/Processed objects/plant_list_dual.rds") # Original code using the whole raw data
saveRDS(df$resource %>% unique, "Data/Processed objects/subset/plant_list_dual.rds")


# Save visitor abundances
insect_ab<-
  polli.div_merged %>% 
  lapply(., function(x){
    lapply(x, function(y){ y %>%
        group_by(from) %>%
        summarize(insect_abundance=sum(weight))
      
    })})


# Save sampling intensity data per site and period
#saveRDS(samp_intensity3, "Data/Processed objects/samp_intensity3_dual.rds") # Original code using the whole raw data
saveRDS(samp_intensity3, "Data/Processed objects/subset/samp_intensity3_dual.rds")


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
  iwint$Int_type2=="mutualistic"] )
### antagonists
table(iwint$con_order[
  iwint$Int_type2=="antagonistic"] )


iwint %>% dplyr::select(year, site_code, record_period, Int_type2, consumer)->smdf
smdf<-smdf%>% unique

smdf %>% group_by(year, site_code, record_period) %>%
  summarise(mutualistic = sum(Int_type2=="mutualistic"),
            antagonistic = sum(Int_type2=="antagonistic"))->smdf

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

#saveRDS(full_data, "Data/Processed objects/degree_unweighted_subset_dual.rds") # Original code using the whole raw data
saveRDS(full_data, "Data/Processed objects/subset/degree_unweighted_dual.rds") 


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
}) #%>% bind_rows(., .id="name")


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

#saveRDS(full_data2, "Data/Processed objects/degree_weighted_dual.rds") # Original code using the whole raw data
saveRDS(full_data2, "Data/Processed objects/subset/degree_weighted.rds") 

#------------------------------------------------------------------------------#
#                                                                              #
#              10. Generating null models for modularity analysis              #
#                                                                              #
#------------------------------------------------------------------------------#

# This part may take some minutes

AM_mat<-polli.div_merged%>%
  map_depth(2, edge_list_to_bipartite)

str(AM_mat, max.level=1)

AM_mat$`2014.A.1.morning`

# Creating a function to generate null models based on Patefield's algorithm (r2d method in bipartite package)
create_null_models <- function(N = 1000) {
  AM_mat_null <- AM_mat %>%
    map_depth(2, ~ {
      tryCatch({
        
        null_matrices <- bipartite::nullmodel(.x*2, method = "r2d", N = N)
        
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
plan(multisession, workers = 8)

with_progress({
  p <- progressor(along = null_layers)
  
  null_multi <- future_imap(null_layers, function(layer_set, u) {
    p()  # Update progress
    map(1:length(layer_set$antagonistic), function(i) {
      lay_pairs <- list(
        antagonistic = layer_set$antagonistic[[i]], 
        mutualistic = layer_set$mutualistic[[i]]
      )
      
      create_multilayer_network(
        list_of_layers = lay_pairs,
        layer_attributes = layer_attrib[[u]], 
        interlayer_links = interlayer_full[[u]],
        bipartite = TRUE,
        directed = FALSE
      )
    })
  }, .options = furrr_options(seed = TRUE))
})

plan(sequential)

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


#saveRDS(map_null,"Data/Processed objects/map_null_dual.rds") # Original code using the whole raw data
saveRDS(map_null,"Data/Processed objects/subset/map_null_dual.rds") 


#------------------------------------------------------------------------------#
#                                                                              #
#                         11. Calculate observed modularity                    #
#                                                                              #
#------------------------------------------------------------------------------#

set.seed(24)
map_all<-lapply(multi.polli,run_infomap_multilayer, flow_model="undirected")
#saveRDS(map_all,"Data/Processed objects/map_all_dual.rds") # Original code using the whole raw data
saveRDS(map_all,"Data/Processed objects/subset/map_all_dual.rds") 


#------------------------------------------------------------------------------#
#                                                                              #
#         12. Preparing data for stability and perturbation analysis           #
#                                                                              #
#------------------------------------------------------------------------------#

polli$period<-sub("^[0-9]+\\.", "", polli$record_period)
polli$network_ID<-paste0(polli$year, ".", polli$site_code)

polli_red<-polli[-which(polli$year==2015|polli$year==2017),]

p_abun<-df2_clean %>% bind_rows(., .id="network_ID") 
p_abun<-p_abun[-which(str_detect(p_abun$network_ID, "2015")|str_detect(p_abun$network_ID, "2017")),]

#write.csv(p_abun, "Data/Processed objects/plant_abundance_F.csv", row.names = FALSE) # Original code using the whole raw data
write.csv(p_abun, "Data/Processed objects/subset/plant_abundance_F.csv", row.names = FALSE) 
#write.csv(polli_red, "Data/Processed objects/polli_red_F.csv", row.names = FALSE) # Original code using the whole raw data
write.csv(polli_red, "Data/Processed objects/subset/polli_red_F.csv", row.names = FALSE)