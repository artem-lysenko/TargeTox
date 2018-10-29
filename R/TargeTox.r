####################################################################################
#
# Method can be run with a call to TargeTox function which takes the following arguments:
# (1) a vector of String ids of protein targets 
# (2-4) route of administration values (administration.oral, administration.parenteral, administration.topical) where
#  values are 0 - not administered via that route; 1 - administered and 2 - unclassified; all defined in equivalent way to ChEMBL database.
# (5-6) Lower and upper plasma protein binding values in a 0-100 range (can be missing)
#
# The method returns a single TargeTox score, with higher score indicating higher toxicity risk. Note that these
# returned values are not probabilities and will not necessarily be within 0-1 range. The optimal cut-offs during
# development was in the region of -1.3 to -1.5 (experimentally confirmed targets from ChEMBL and DrugBank), though we recommend
# that users perform their own calibration if using definitions of targets not consistent with ours (e.g. if using
# computationally predicted drug target sets with high false-positive rates).
#
####################################################################################
library(catboost)
library(sets)


load(paste0(dirname(sys.frame(1)$ofile),"/resources.RData"))
mod = catboost.load_model(paste0(dirname(sys.frame(1)$ofile),"/TargeTox.cbm"))

getNonredundat_v2 = function(input){
    output = as.set(input)
    for(term in input){
        output = set_complement(ancMap[[term]], output)    
    }
    return (as.character(output))
}

functionalImpact = function(targets){
    return(sum(ic_lookup[getNonredundat_v2(string2go[which(string2go$string %in% targets),'go'])]))
}

minDsd = function(targets){
    return(sapply(1:ncol(dsdMatrix), FUN=function(ce){
        min(dsdMatrix[which(rownames(dsdMatrix) %in% targets) , ce])
    }))    
}


TargetTox = function(targets, administration.oral,
                     administration.parenteral, administration.topical,
                     protein.binding.lower=NA, protein.binding.upper=NA){
    if(is.na(administration.oral) & is.na(administration.parenteral) & is.na(administration.topical)){
        stop("At least one route of administration should be set.")    
    }
    input = c(minDsd(targets), functionalImpact(targets), administration.oral,
        administration.parenteral, administration.topical,
        protein.binding.lower, protein.binding.upper)

    input = as.data.frame(t(as.matrix(input)))

    #colnames(input) = c('9606.ENSP00000438284', '9606.ENSP00000232975', '9606.ENSP00000366416',
    #    '9606.ENSP00000315654', '9606.ENSP00000261712', '9606.ENSP00000308549', '9606.ENSP00000326042',
    #    '9606.ENSP00000281821', '9606.ENSP00000226021', '9606.ENSP00000328625', '9606.ENSP00000447149',
    #    '9606.ENSP00000335592', 'drug_ic', 'oral', 'parenteral', 'topical', 'lower', 'higher') 

	colnames(input) = c('n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11',
        'n12', 'drug_ic', 'oral', 'parenteral', 'topical', 'lower', 'higher') 
	
    input$oral = factor(input$oral, levels = c('0', '1', '2')) 
    input$parenteral = factor(input$parenteral, levels = c('0', '1', '2')) 
    input$topical = factor(input$topical, levels = c('0', '1', '2')) 

    types = c('numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric',
              'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'factor', 'factor', 'factor',
              'numeric', 'numeric')
    lp1 = catboost.load_pool(input, column_description = types)

    result = catboost.predict(mod, lp1)
    return(result)
}

# Test example
#ids = c('9606.ENSP00000338072', '9606.ENSP00000237837')
#administration.oral = 0
#administration.parenteral = 1
#administration.topical = 0
#protein.binding.lower = NA
#protein.binding.upper = NA

#TargetTox(ids, administration.oral, administration.parenteral, administration.topical,
#                     protein.binding.lower, protein.binding.upper)


