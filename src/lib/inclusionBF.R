inclusionBF <- function(BFobj, effects = NULL, models = "all") {
  # Computes the "inclusion Bayes factor" for a particular effect from the Bayes factors for individual models.
  #
  # Args:
  #   BFobj: an object of class 'BFBayesFactor', which is the default output of functions in the BayesFactor package
  #   effect: a character vector with the name(s) of the effect(s) you want to obtain the inclusion Bayes factor of.
  #           for example, "age", or c("age", "gender:length").
  #           Omitting this input will return the inclusion Bayes factors for all effects.
  #   models: either "all" to compute the standard inclusion Bayes factor [default], or
  #           "matched" to compute across matched models (see below).
  #
  # Returns:
  #   a data frame with the inclusion Bayes factors for all / the specified effects
  #
  # The default inclusion Bayes factor (option "all") is a concept from the JASP stats package: https://jasp-stats.org/
  # Briefly, it compares all models that include the effect of interest vs. all models that do not.
  # For examples and a conceptual explanation, see example 5 in: 
  #   "Bayesian inference for psychology. Part II: Example applications with JASP" https://doi.org/10.3758/s13423-017-1323-7
  #
  # The Bayes factor across matched models (option "matched") is currently being implemented in JASP.
  # It was conceptualized by Sebastiaan Mathot, see:
  #   https://www.cogsci.nl/blog/interpreting-bayesian-repeated-measures-in-jasp (and associated posts on the forum)
  # It is more selective in the set of models that is being compared than the standard inclusion BF.
  # Briefly, the inclusion BF across matched models compares:
  #   - all models that include the effect of interest, but NO interactions with the effect of interest, VERSUS
  #   - the models that result from stripping the effect of interest from this set of models
  
  # L.C. Reteig, 03/08/2017
  # Based on code by Florian Sense: http://pastebin.com/Nfmw4JDs
  
  BFs <- extractBF(BFobj)[,1:2] # convert to data frame, keep only first two columns
  
  if (is.null(effects)) { # if no effect was specified
    effects <- gsub(" ","",strsplit(rownames(BFs)[which.max(nchar(rownames(BFs)))], "\\+")[[1]]) # get the names of all the effects in the model
    
    if ((unlist(attributes(BFobj))$denominator@dataTypes == "random")) { # if there's a random factor
      randName <- names(unlist(attributes(BFobj))$denominator@dataTypes == "random") # get its name
      effects <- effects[!(effects %in% randName)] # and strip it from the list of effects
    }
  }
  
  BFs <- rbind(BFs,data.frame(bf = 1, error = 0, row.names = "null")) # add the null model, which by definition has a BF10 of 1
  
  inclBF <- data.frame(effect = effects, Bayes.factor = NA) # initialize data frame
  for (term in effects) {
    
    idx <- grepl(term, rownames(BFs)) # search for rows containing the effect of interest
    
    if (all(!idx)) stop(paste("The effect", term, "does not appear in any model. Inclusion BF cannot be calculated."))
    
    if (models == "all") {
      
      # Calculate prior inclusion probability:
      # assuming all models are equally likely, this is simply the proportion of models that include the effect of interest
      prior <- sum(idx) * (1/nrow(BFs))
      
      # Calculate posterior inclusion probability: 
      # this is the sum of the posterior probability of the models that include the effect of interest.
      # BayesFactor does not directly output posterior probabilities, but they can be calculated from the Bayes Factors,
      # because the posterior probabilities must sum to 1, and are proportional to the BF10
      post <- sum(BFs$bf[idx]) / sum(BFs$bf)
      
      if (post == 1) warning("The posterior inclusion probability is so high that R rounds it to 1. The inclusion BF will be Inf.")
      
      # Convert to odds
      priorOdds <- (prior / (1 - prior)) # ratio of the posterior inclusion probability versus the non-inclusion probability
      postOdds <- (post / (1 - post)) # ratio of the prior inclusion probability versus the non-inclusion probability
      
      # The inclusion Bayes factor is the posterior odds divided by the prior odds (i.e. the change from prior to posterior inclusion odds)
      inclBF$Bayes.factor[inclBF$effect %in% term] <- postOdds / priorOdds
      
    } else if (models == "matched") {
      
      # With-models:
      # If the effect of interest is a main effect:
      # - contain the term itself, so for example "gender"
      # - do not have interactions with the effect of interest, so do not contain terms like age:gender
      # If the effect of interest is an interaction, the same rules apply:
      # - contain the term itself, so for example "gender:length"
      # - do not have interactions with the effect of interest, so do not contain terms like "age:gender:length"
      
      # Calculate all possible permutations of higher-order interactions
      effect_factors <- strsplit(term,":")[[1]] # names of all main effects
      effect_factors <- c(effect_factors, "\\S+") # add a "non-whitespace" regexp
      perms <- matrix(effect_factors[permutations(length(effect_factors))],ncol = length(effect_factors))
      
      int_mat <- matrix(paste0(perms,":"),nrow = nrow(perms)) # add a: to each effect
      # append a "|" to each permutation; take it away for the last permutation
      int_mat[,ncol(perms)] <- matrix(sub(":","|",int_mat[,ncol(perms)]),nrow = nrow(perms))
      int_mat[nrow(int_mat),ncol(int_mat)] <- sub("\\|","",int_mat[nrow(int_mat),ncol(int_mat)])
      # collapse everything into a single string to make one regular expression
      reg_exp_int <- paste(t(int_mat),collapse = "") # regexp for all possible higher-order interactions
      
      # With models:
      withModelIdx <- grepl(term, rownames(BFs)) & # must have term itself AND
        !(grepl(reg_exp_int, rownames(BFs))) # no interactions with the term
      
      # Without-models:
      # The set of without-models is the set of with-models, with the term of interest stripped from each model
      
      withModels <- rownames(BFs[withModelIdx,]) # names of all the with-models
      sS <- paste0(term, " \\+ ", "|", " \\+ ", term) # regular expression for removing the term of interest, e.g. "+ age" OR "age +"
      withoutModels = sub(sS,"",withModels) # remove these terms, to yield the set of without-models
      
      # Add the null model to the set of without-models if necessary
      if ( !(grepl(":", term)) ) { # if the effect of interest is a main effect
        
        if (unlist(attributes(BFobj))$denominator@dataTypes == "fixed") { # if all factors are fixed
          # for the with-model that contains only the main effect of interest, the "null" model is its without-model
          withoutModels[withModels %in% term] <- "null"
          
        } else if ((unlist(attributes(BFobj))$denominator@dataTypes == "random")) { # if there's a random factor in the denominator
          randName <- names(unlist(attributes(BFobj))$denominator@dataTypes == "random") # get its name
          # for the with-model that contains only the main effect of interest PLUS the random factor, the "null" is its without-model
          withoutModels[withModels %in% paste(term, "+", randName)] <- "null"
        }
      }
      
      withoutModelIdx <- rownames(BFs) %in% withoutModels # indices of without-models
      
      postWith <- sum(BFs$bf[withModelIdx]) / sum(BFs$bf) # posterior inclusion probability for with-models
      postWithout <- sum(BFs$bf[withoutModelIdx]) / sum(BFs$bf) # posterior inclusion probability for without-models
      if (postWithout == 1) warning("The posterior inclusion probability for without-models is so high that R rounds it to 1. The inclusion BF will be Inf.")
      # Note that we don't need priors or odds ratios, as there are always as many with-models as without-models
      
      # The inclusion Bayes Factor across matched models is the ratio between with- and without-model posteriors
      inclBF$Bayes.factor[inclBF$effect %in% term] <- postWith / postWithout 
    }
    
  }
  
  inclBF
}

# Function for generating permutations, from
# https://stackoverflow.com/a/20199902
permutations <- function(n){
  if (n == 1) {
    return(matrix(1))
  } else {
    sp <- permutations(n - 1)
    p <- nrow(sp)
    A <- matrix(nrow = n*p,ncol = n)
    for (i in 1:n) {
      A[(i - 1)*p + 1:p,] <- cbind(i,sp + (sp >= i))
    }
    return(A)
  }
}