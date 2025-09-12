# Custom Function for Metabolomics Data Analysis 
# Ghiwa Makke
#
#

# Function to add spaces between chemical elements and quantities, but keep two-letter elements intact
add_spaces_to_formula <- function(formula) {
  # Add spaces between numbers and the next element
  spaced_formula <- gsub("([0-9])([A-Z])", "\\1 \\2", formula)
  # Add spaces between elements but avoid separating two-letter elements
  spaced_formula <- gsub("([A-Z][a-z]*)(?=[A-Z])", "\\1 ", spaced_formula, perl = TRUE)
  return(spaced_formula)
} 

# Function by Christian Ayala-Ortiz-------------------------------------------------------------------------
get_elements <- function(df){
  ## Obtain element list in the formula
  element_list <- gsub('[[:digit:]]+', '', df$Formula)
  element_list <- paste(element_list, collapse = ' ')
  element_list <- str_split(element_list, ' ')[[1]]
  element_list <- unique(element_list)
  
  return(element_list)
}
# Function by Christian Ayala-Ortiz-------------------------------------------------------------------------------
separate_formula <- function(df){
  # This function will split the Formula column
  # into columns with the number of each element
  
  ## Get the elements that make each of the compounds
  element_list <- get_elements(df)
  
  ## Get formula in a df were results will be stored
  result_df <- select(df, Formula)
  
  ## Separate Formula into the elements
  new_df <- separate(result_df, 
                     Formula,
                     into = {{element_list}},
                     sep = ' ')
  
  ## Initialize a temporal df and an accumulator for the for loops
  temp_df <- tibble(.rows = nrow(new_df))
  j <- 1
  
  ## For loop to obtain the coefficients of each element
  for(el in element_list){
    my_exp <- paste0(el, "\\d\\d|", el, "\\d")
    for(i in 1:length(element_list)){
      temp_df[,i] <- ifelse(grepl(el, deframe(new_df[,i])), # because new_df is a tibble, deframe allows it to be sliced as a vector
                            ifelse(grepl(my_exp, deframe(new_df[,i])), str_replace(deframe(new_df[,i]), paste0(el), ''), 1), 0)
      temp_df[,i] <- as.numeric(unlist(temp_df[,i])) # unlist produces atomic components that can be changed into numeric values
    }
    temp_df$sum <- rowSums(temp_df)
    result_df[,j + 1] <- temp_df$sum
    j <- j + 1
    temp_df <- tibble(.rows = nrow(new_df))
  }
  
  ## Put back the names of the elements in each column
  colnames(result_df) <- c('Formula', element_list)
  
  ## Merge with original matrix
  df <- left_join(df, result_df, by = "Formula")
  
  df <- distinct(df)
  
  return(df)
  
}

# Function by Christian Ayala-Ortiz-------------------------------------------------------------------------
calc_ratios_n_idxs <- function(df){
  # This function will calculate H/c and O/C ratios
  # as well as other thermodynamics index
  
  df$C <- as.numeric(df$C)
  df$H <- as.numeric(df$H)
  df$O <- as.numeric(df$O)
  df$N <- as.numeric(df$N)
  df$P <- as.numeric(df$P)
  df$S <- as.numeric(df$S)
  
  ## Get ratios
  df <- df %>% 
    mutate(H_to_C = H / C) %>% 
    mutate(O_to_C = O / C)
  
  ## Calculate thermodynamic indices
  df <- df %>% 
    mutate(NOSC = -((4*C + H - 3*N - 2* O + 5*P - 2*S) / C) + 4) %>% 
    mutate(GFE = 60.3 - 28.5*NOSC) %>% 
    mutate(DBE = 1 + 0.5 * (2*C - H + N + P)) %>% 
    mutate(DBE_O = DBE - O) %>% 
    mutate(AI = (1 + C - O - S - ((H + P + N) * 0.5)) / (C - O - S - N - P)) %>% 
    mutate(AI_mod = (1 + C - (O * 0.5) - S - ((H + P + N) * 0.5)) / (C - (O * 0.5) - S - N - P)) %>% 
    mutate(DBE_AI = 1 + C -O -S - (0.5 * (H + N + P)))
}

# Function by Christian Ayala-Ortiz-------------------------------------------------------------------------
calc_classes <- function(df){
  df <- df %>% 
    mutate(Class = ifelse(between(O_to_C, 0, 0.3)&between(H_to_C, 1.5, 2.5), 'Lipid',
                          ifelse(between(O_to_C, 0, 0.125)&between(H_to_C, 0.8, 1.5), 'Unsaturated HC',
                                 ifelse(between(O_to_C, 0, 0.95)&between(H_to_C, 0.2, 0.8), 'Condensed HC',
                                        ifelse(between(O_to_C, 0.3, 0.55)&between(H_to_C, 1.5, 2.3), 'Protein',
                                               ifelse(between(O_to_C, 0.55, 0.7)&between(H_to_C, 1.5, 2.2), 'Amino Sugar',
                                                      ifelse(between(O_to_C, 0.7, 1.5)&between(H_to_C, 1.5, 2.5), 'Carbohydrate',
                                                             ifelse(between(O_to_C, 0.125, 0.65)&between(H_to_C, 0.8, 1.5), 'Lignin',
                                                                    ifelse(between(O_to_C, 0.65, 1.1)&between(H_to_C, 0.8, 1.5), 'Tannin', 'Other')))))))))
  return(df)
}


#Normalization functions ---------------
#Pareto Scaling:
PS_helper <- function(x) {
  (x - mean(x)) / sqrt(sd(x, na.rm = T))
}	

pareto_scale <- function(x){
  mtb_scaled <- data.frame(apply(x, 2, PS_helper))
  return(mtb_scaled)
}

#Auto Scaling:
AS_helper <- function(x) {
  (x - mean(x)) / sd(x, na.rm = T)
} 

auto_scale <- function(x){
  mtb_scaled <- apply(x, 2, AS_helper) 
  return(mtb_scaled)
}

#Log Transformation Functions:
log_helper <- function(x, min.val) {
  log2((x + sqrt(x ^ 2 + min.val ^ 2)) / 2)
}

#Log Scaling:
log_transform <- function(x){
  x_nz <- x[ ,which(apply(x, 2, sum) != 0)]
  min.val <- min(abs(x_nz[x_nz!=0]))/10
  x_log_trans <- data.frame(apply(x_nz, 2, log_helper, min.val))
  return(x_log_trans)
}

calc_log2fc <- function(x, factor, groups){
  if(length(groups) != 2){
    stop('Must Be Exactly 2 Groups')
  }
  group1 <- x[x[factor] == groups[1], ][['area']]
  group2 <- x[x[factor] == groups[2], ][['area']]
  l2fc <- log2(mean(group2)/mean(group1))
  return(l2fc)
}
