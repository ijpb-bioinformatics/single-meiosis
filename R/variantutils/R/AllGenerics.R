#' @title ConvertToVcf Method
#' @description Convert VCFlike data to VCF 
#' @details This function is called by the class VCFlikeParser to set the slot VCF of an object of class VCFlike from
#' the Df and MetaDf slots of this object.
#' @param VCFlike A VCFlike object
#' @return VCF
#' @name VCFlike
#' @rdname VCFlike-class

setGeneric(
    name="ConvertToVcf",
    def = function( object, ... ){standardGeneric("ConvertToVcf")}
    )

#' @title Method Genotype
#' @description Annotate a VCF by adding genotype information to the VCF object of a VCFlike.
#' @return Add or replace the geno slot of a VCF object 
#' @name SM-class
#' @rdname SM-class

setGeneric(
  name="Genotype",
  def = function( object, ... ){standardGeneric("Genotype")}
)

#' @title Method Filter
#' @description Annotate a VCF by adding filtering information to the VCF object of a VCFlike.
#' @return Add or replace the filter slot of a VCF object 
#' @name SM-class
#' @rdname SM-class

setGeneric(
  name="Filter",
  def = function( object, ... ){standardGeneric("Filter")}
)


setGeneric("AnnotateVCF", 
          def = function(object, ...)
             standardGeneric("AnnotateVCF")
)

#' @title Method FindNCO
#' @description Find putativ NCO events.
#' @return Add or replace the DfOfNCO slot of an object of Class Tetrad. 
#' @name Tetrad-class
#' @rdname Tetrad-class

setGeneric(name="FindNCO",      
          def = function(object){standardGeneric("FindNCO")})


#' @title Method FindCO
#' @description Find putativ CO events.
#' @return Add or replace the DfOfCO slot of an object of Class Tetrad. 
#' @name Tetrad-class
#' @rdname Tetrad-class

setGeneric(name="FindCO",      
          def = function(object, ...){standardGeneric("FindCO")})


#' @title Method setDfOfNCO
#' @return Add or replace the DfOfNCO slot of a Tetrad object 
#' @exportMethod setDfOfNCO
#' @name Tetrad-class
#' @rdname Tetrad-class

setGeneric(name="setDfOfNCO<-",      
           def = function(object,...){standardGeneric("setDfOfNCO<-")})


#' @title Method setDfOfCO
#' @return Add or replace the DfOfNCO slot of a Tetrad object 
#' @exportMethod setDfOfCO
#' @name Tetrad-class
#' @rdname Tetrad-class

setGeneric(name="setDfOfCO<-",      
           def = function(object,...){standardGeneric("setDfOfCO<-")})


#' @title Method writeTetrad
#' @description Export the Tetrad in VCF format and export the two data frame dfOfCO and
#' dfOfNCO into a text file.
#' @name SM-class
#' @rdname SM-class

setGeneric(name="writeTetrad",      
          def = function(object, ...){standardGeneric("writeTetrad")})


#setGeneric(
#    name="countMismatches",
#    def = function( object, tags){standardGeneric("countMismatches")}
#    )

setGeneric(
	name = "countMismatches",
	def = function(object, tags, ...) {
		df <- standardGeneric("countMismatches")
		if(!(is(df, "data.frame"))) 
			stop("countMismatches methods must return data frame objects")
		df
	})

setGeneric(
    name = "countPrimaryTag",
    def = function(object, tags, ...) {
        df <- standardGeneric("countPrimaryTag")
        if(!(is(df, "data.frame")))
            stop("countPrimaryTag method must return data frame objects")
        df
    })

setGeneric(
		   name = "setCounts<-",
		   def = function(object,value){standardGeneric("setCounts<-")})

