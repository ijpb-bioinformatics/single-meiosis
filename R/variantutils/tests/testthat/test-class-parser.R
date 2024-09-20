library(testthat)
library(VariantAnnotation)
library(stringr)
library(seqinr)

TAIR10 <- read.fasta("/data/SEQUENCES/GENOME/Ath-Col0-tair10-WG/Ath-Col0-tair10-WG.mfa")
RefGenomes <- list("TAIR10"=TAIR10)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir("../../R/")

test_that("TestInstanceParserClassRefPol",{
  testParser <- new("VCFlikeParser",Format="RefPol",File="test-class-SMsdi.txt",RefGenome="TAIR10",Sample="test")
  writeVcf(testParser@VCFlike@VCF,file="test-class-RefPol.vcf")
  expect_that(class(testParser),matches("VCFlikeParser"))
})

AnnotateVCF(testParser@VCFlike@VCF, vcfField="IMR_DENOM_SCORE", vcfFieldType=c("info"), FilterType=c("chararacter"), 
FilterName="IMR_DENOM_SCORE", FilterValue="2", FilterDescription="score",  append=TRUE)

test_that("TestInstanceParserClassSM",{
  testParser <- new("VCFlikeParser",Format="SM",File="../../data/test.SM",RefGenome="coller-0",Sample="VersaillesTetrad4IndividuM1")
  expect_that(class(testParser),matches("VCFlikeParser"))
})

test_that("TestInstanceParserClassSM",{
  testParser <- new("VCFlikeParser",Format="SM",File="../../data/test.SM",RefGenome="coller-0",Sample="VersaillesTetrad4IndividuM1")
  expect_that(class(testParser),matches("VCFlikeParser"))
})


# One Substitution report test

test_that("SubstitutionLigne3DansTestSdi9",{
    testParser <- new("VCFlikeParser",Format="sdi-9",File="../../data/test.sdi9",RefGenome="TAIR10",Sample="oy-0")
    #start
    expect_that(start(rowData(testParser@VCFlike@VCF))[3],equals(107))
    #end
    expect_that(end(rowData(testParser@VCFlike@VCF))[3],equals(107))
    #width
    expect_that(width(rowData(testParser@VCFlike@VCF))[3],equals(1))
    #reference_base
    expect_that(as.vector(ref(testParser@VCFlike@VCF)[3]),matches("C"))
    #alternate_base
     expect_that(as.vector(alt(testParser@VCFlike@VCF)[[3]]),matches("T"))
    })

# java -jar /usr/local/src/snpEff_development/snpEff.jar eff -c /usr/local/src/snpEff_development/snpEff.config -i vcf -o vcf athalianaTair10 testsdi9export.vcf
# One Deletion report test

test_that("DeletionLigne17DansTestSdi9",{
    testParser <- new("VCFlikeParser",Format="sdi-9",File="../../data/test.sdi9",RefGenome="TAIR10",Sample="oy-0")
    #start
    expect_that(start(rowData(testParser@VCFlike@VCF))[17],equals(701))
    #end
    expect_that(end(rowData(testParser@VCFlike@VCF))[17],equals(700))
    #width
    expect_that(width(rowData(testParser@VCFlike@VCF))[17],equals(0))
    #reference_base
    expect_that(as.vector(ref(testParser@VCFlike@VCF)[17]),matches("T"))
    #alternate_base
     expect_that(as.vector(alt(testParser@VCFlike@VCF)[[17]]),matches(""))
    })


test_that("DeletionLigne21DansTestSdi9",{
    testParser <- new("VCFlikeParser",Format="sdi-9",File="../../data/test.sdi9",RefGenome="TAIR10",Sample="oy-0")
    #start
    expect_that(start(rowData(testParser@VCFlike@VCF))[21],equals(800))
    #end
    expect_that(end(rowData(testParser@VCFlike@VCF))[21],equals(799))
    #width
    expect_that(width(rowData(testParser@VCFlike@VCF))[21],equals(0))
    #reference_base
    expect_that(as.vector(ref(testParser@VCFlike@VCF)[21]),matches("CTTT"))
    #alternate_base
     expect_that(as.vector(alt(testParser@VCFlike@VCF)[[21]]),matches(""))
    })


# One Insertion report test

test_that("InsertionLigne2DansTestSdi9",{
    testParser <- new("VCFlikeParser",Format="sdi-9",File="../../data/test.sdi9",RefGenome="TAIR10",Sample="oy-0")
    #start
    expect_that(start(rowData(testParser@VCFlike@VCF))[2],equals(35))
    #end
    expect_that(end(rowData(testParser@VCFlike@VCF))[2],equals(42))
    #width
    expect_that(width(rowData(testParser@VCFlike@VCF))[2],equals(8))
    #reference_base
    expect_that(as.vector(ref(testParser@VCFlike@VCF)[2]),matches(""))
    #alternate_base
     expect_that(as.vector(alt(testParser@VCFlike@VCF)[[2]]),matches("GATCCTT"))
    })


test_that("TestMethodGenotypeSM",{
  testParser <- new("VCFlikeParser",Format="SM",File="../../data/test.SM",RefGenome="coller-0",Sample="M1")
  
  testParser@VCFlike <- Genotype(testParser@VCFlike,0.003)
  writeVcf(testParser@VCFlike@VCF,file="testSM.txt")
})

test_that("TestMethodFilterSM",{
  testParser <- new("VCFlikeParser",Format="SM",File="../../data/test.SM",RefGenome="coller-0",Sample="VersaillesTetrad4IndividuM1")
  
  testParser@VCFlike@VCF <- AnnotateVCF(testParser@VCFlike@VCF,vcfField="COL_DP",vcfFieldType="info", FilterType="numeric",FilterName="COL_DP",FilterValue=20,FilterSign="more",FilterDescription="Depth on Columbia Genome", append=TRUE)
  
  writeVcf(testParser@VCFlike@VCF,file="testSM_Annotate.txt")
})

 test_that("TestClassTetrad",{
   testParser1 <- new("VCFlikeParser",Format="SM",File="../../data/test.SM",RefGenome="coller-0",Sample="M1")
   testParser2 <- new("VCFlikeParser",Format="SM",File="../../data/test2.SM",RefGenome="coller-0",Sample="M2")
   testParser3 <- new("VCFlikeParser",Format="SM",File="../../data/test.SM",RefGenome="coller-0",Sample="M3")
   testParser4 <- new("VCFlikeParser",Format="SM",File="../../data/test2.SM",RefGenome="coller-0",Sample="M4")
   
   VCFL <- list(testParser1,testParser2,testParser3,testParser4)
   LVCF <- lapply(VCFL, function(x){
     x@VCFlike <- Genotype(x@VCFlike,0.003)
     return(x)
   })
   VCFTetrad <- new("Tetrad",LVCF)  
 })

test_that("TestFindNCO",{
  testParser1 <- new("VCFlikeParser",Format="SM",File="../../data/test3.SM",RefGenome="coller-0",Sample="M1")
  testParser2 <- new("VCFlikeParser",Format="SM",File="../../data/test2.SM",RefGenome="coller-0",Sample="M2")
  testParser3 <- new("VCFlikeParser",Format="SM",File="../../data/test2.SM",RefGenome="coller-0",Sample="M3")
  testParser4 <- new("VCFlikeParser",Format="SM",File="../../data/test2.SM",RefGenome="coller-0",Sample="M4")
  
  VCFL <- list(testParser1,testParser2,testParser3,testParser4)
  LVCF <- lapply(VCFL, function(x){
    x@VCFlike <- Genotype(x@VCFlike,0.003)
    return(x)
  })
  VCFTetrad <- new("Tetrad",LVCF)  
  VCFTetrad <- FindNCO(VCFTetrad)
})

test_that("TestPlot",{
  testParser1 <- new("VCFlikeParser",Format="SM",File="../../data/test3.SM",RefGenome="coller-0",Sample="M1")
  testParser2 <- new("VCFlikeParser",Format="SM",File="../../data/test2.SM",RefGenome="coller-0",Sample="M2")
  testParser3 <- new("VCFlikeParser",Format="SM",File="../../data/test2.SM",RefGenome="coller-0",Sample="M3")
  testParser4 <- new("VCFlikeParser",Format="SM",File="../../data/test2.SM",RefGenome="coller-0",Sample="M4")
  
  VCFL <- list(testParser1,testParser2,testParser3,testParser4)
  LVCF <- lapply(VCFL, function(x){
    x@VCFlike <- Genotype(x@VCFlike,0.003)
    return(x)
  })
  VCFTetrad <- new("Tetrad",LVCF)  
  plot(VCFTetrad)
})


test_that("TestInstanceParserClassJGI-aa",{
  testParserJGIaa <- new("VCFlikeParser",Format="JGI-aa",File="JGI_test_aa.txt",RefGenome="TAIR10",Sample="Ac-0")
  writeVcf(testParserJGIaa@VCFlike@VCF,file="JGIaa_test_aa_tmp.vcf")
  expect_that(class(testParserJGIaa),matches("VCFlikeParser"))
})

test_that("TestInstanceParserClassSDI-9",{
  testParserSdi9 <- new("VCFlikeParser",Format="sdi-9",File="sdi-9_test_Bur0.txt",RefGenome="TAIR10",Sample="bur-0")
  writeVcf(testParserSdi9@VCFlike@VCF,file="sdi-9_test2_Bur0.vcf")
  expect_that(class(testParserSdi9),matches("VCFlikeParser"))
})


test_that("TestInstanceParserClassJGI-indel",{
  testParserJGIindel <- new("VCFlikeParser",Format="JGI-indel",File="JGI_test_indel.txt",RefGenome="TAIR10",Sample="Ac-0")
  writeVcf(testParserJGIindel@VCFlike@VCF,file="JGI_test_indel.vcf")
  expect_that(class(testParserJGIindel),matches("VCFlikeParser"))
})

test_that("TestInstanceParserVAST",{
  testParserVAST <- new("VCFlikeParser",Format="VAST",File="VAST_test_Cvi.txt",RefGenome="TAIR10",Sample="Cvi")
  writeVcf(testParserVAST@VCFlike@VCF,file="VAST_test_Cvi.vcf")
  expect_that(class(testParserVAST),matches("VCFlikeParser"))
})

test_that("TestInstanceParserSALK",{
  testParserSALK <- new("VCFlikeParser",Format="SALK",File="SALK_test_etna-2.txt",RefGenome="TAIR10",Sample="etna-2")
  writeVcf(testParserSALK@VCFlike@VCF,file="SALK_test_etna-2.vcf")
  expect_that(class(testParserSALK),matches("VCFlikeParser"))
})



