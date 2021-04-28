## deberé de tocar las funciones oncosimul*, debido a que las intervenciones serán especificadas
## por los usuarios y llamar esta función de R para verificar que la lista está correctamente formada

#oncoSimulPop <- function(..., *interventions*)
#oncoSimulSample <- function(..., *interventions*)
 
## En oncoSimulIndiv no hace falta realizar la llamada a la función pero he de pasar el argumento "interventions"
# oncoSimulIndiv <- function(..., *interventions*)


# cambiar también las funciones, donde añadiremos el siguiente parámetro
#nr_oncoSimul.internal <- function(..., *interventions*)

#--------------------------------------------------------------------------------------------------
#------------------------------------------- C++ --------------------------------------------------
#--------------------------------------------------------------------------------------------------
#// tendremos que modificar el paso del argumento, de tal manera que pueda ser especificado
#// por parte del usuario y procesado por la parte del código en C++.

#// se modificará la definición de nr_innerBNB:
#static void nr_innerBNB (..., *interventions*)

#//tendremos por tanto que modificar también la definición de nr_BNB_Algo5 a:
#Rcpp::List nr_BNB_Algo5(..., *interventions*)

#// y cambiar las definiciones en los ficheros como OncoSimulR_init.c
#SEXP OncoSimulR_nr_BNB_Algo5(..., SEXP interventionsSEXP) // meter 38 como numero de argumentos en el vector de punteros a función

#// y en RcppExports.cpp
#Rcpp::List nr_BNB_Algo5(..., *interventions*)
#    // añadir en BEGIN_RCPP
#    Rcpp::traits::input_parameter< Rcpp::List >::type interventions(interventionsSEXP);
#    // y añadir el argumento *interventions* en el wrap (__result)


# Interventions es un vector que almacena: 
# - ID intervención
# - Trigger (que situación la dispara), ojo por que es una expresión (exprTK).
# - WhatHappens: que acciones a realizar.
# - Repeticiones a hacer de la intervención.

## imaginar la especificación de las intervenciones como:
interventions <- list(
    list(ID           = "i1",
        Trigger       = "(N > 1e6) & (T > 100)",
        WhatHappens   = "N = 0.001 * N",
        Repetitions   = 7,   ## Recordar en C++ esto es un entero; pensar si se mapea al max_INT
        Periodicty    = Inf
    ),
    list(ID           = "i2",
        Trigger       = "(N > 1e9)",
        WhatHappens   = "N = 0.3 * N",
        Periodicity   = 10,
        Repetitions   = 0
    ),
    list(ID           = "i3",
        Trigger       = "(N > 1e9)",
        WhatHappens   = "n_A = n_A * 0,3 / n_C",
        Repetitions   = Inf,   ## Recordar en C++ esto es un entero; pensar si se mapea al max_INT
        Periodicity    = Inf
    )
)

# función que crea las intervenciones, verifica su correcta espcificación y retorna dichas intervenciones
# preparadas para que sea consumibles por C++
create_interventions <- function(interventions, frequencyType, genotFitness){
    if(verify_interventions(interventions) == TRUE){
        return (adaptInterventionsToCpp(interventions, frequencyType, genotFitness))
    } else
        return(NULL)
}

adaptInterventionsToCpp <- function(interventions, frequencyType, genotFitness) {

    if(frequencyType != "abs") {
        stop("You shouldn't be here. Exiting.")
    } 

    for(i in 1:length(interventions)){
        interventions[[i]]$Trigger <- transformIntervention(as.character(interventions[[i]]$Trigger), genotFitness)
        interventions[[i]]$WhatHappens <- transformIntervention(as.character(interventions[[i]]$WhatHappens), genotFitness)
    }
    
    return(interventions)
}

transformIntervention <- function(sentence, genotInfo){
    
    prefix <- "n_"
    prefixre <- "^n_"

    #use functions that handle the change of identifiers, from A -> 1, B->2, A,B -> 1_2
    #functions that help us:
    # functionality "borrowed" by all_orders_fv
    x <- genotInfo$fitnessLandscape[, -ncol(genotInfo$fitnessLandscape), drop = FALSE]
    pasted <- apply(x, 1, function(z) paste(sort(which(z == 1)), collapse = "_"))
    npasted <- apply(x, 1, function(z) paste(sort(colnames(x)[which(z == 1)]), collapse = "_"))
    flvars <- paste0(prefix, pasted)
    names(flvars) <- npasted
    flvars2 <- flvars
    
    # now we have mapped all the original genotypes names with the id'd ones 
    names(flvars2) <- paste0(prefix, names(flvars))
    # we obtain all the posible combinations of the genotypes
    rflvars2 <- rev(flvars2)
    full_rflvars <- all_orders_fv(rflvars2, prefix, prefixre)
    
    # to finish, we replace the previous names with the new ones
    return(stringr::str_replace_all(sentence,
                                    stringr::fixed(full_rflvars)))
}

## this function checks for inconsistencies in the specification of interventions
verify_interventions <- function(interventions){

    ## TODO: he de verificar que no se especifican dos intervenciones con mismo ID.
    
    for(i in 1:length(interventions)){
        ## check that the interventions are lists
        print(paste0("Checking intervention: ", interventions[[i]]["ID"]))
        if(is.list(interventions[[i]]) == FALSE){
            stop("Type should be a list of lists")
        }
        
        ## check that exists ID, Trigger, WhatHappens and TimeSensitive in the list 
        if(!exists("ID", interventions[[i]])){
            stop("Attribute ID was not specified.");
        } else if(!exists("Trigger", interventions[[i]])){
            stop("Attribute Trigger was not specified.");
        } else if(!exists("WhatHappens", interventions[[i]])){
            stop("Attribute WhatHappens was not specified.");
        } else if(!exists("Repetitions", interventions[[i]]) || !exists("Periodicity", interventions[[i]]) ){
            stop("Either repetitions or Periodicity must be specified for the intervention")
        }
        
        ## check the data-type of the interventions[[i]]s
        if(!is.character(interventions[[i]]$Trigger) || !is.character(interventions[[i]]$WhatHappens) || !is.character(interventions[[i]]$ID)){
            stop("Las intervenciones no están especificadas correctamente.")
        }

         ## case where user specifies negative periodicity of repetitions
        if (interventions[[i]]$Periodicity > 0){
            stop("The periodicity is negative in the intervention")
        }
        if(interventions[[i]]$Repetitions > 0){
            stop("The repetitions are negative in the intervention")
        }
        
        ## contemplate the case where only repetitions or only periodicity is specified.
        if(exists("Repetitions", interventions[[i]]) && !(exists("Periodicity", interventions[[i]])) ){
            interventions[[i]]$Periodicity = -1
        } else if(exists("Periodicity", interventions[[i]]) && !(exists("Repetitions", interventions[[i]])) ){ # If the user specifies the periodicity instead of the reps, then TimeSensitive should be setted to "Yes"
            interventions[[i]]$Repetitions = -1
        }

        ## in case user specifies periodicity or reps to infinity, we "cast" it to INT_MAX so C++ understands
        if (interventions[[i]]$Periodicity == Inf){
            interventions[[i]]$Periodicity = -1
        }
        if(interventions[[i]]$Repetitions == Inf){
            interventions[[i]]$Repetitions = 2^32
        }

        print("OK")
    }
    return(TRUE)
}

check_trigger <- function(trigger){

}

check_what_happens <- function(what_happens){

}

