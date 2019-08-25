########################################################################################
########################################################################################
##############  FUNCTIONS FOR PROCESSING OTU GROUPING USING ############################ 
############## CONTRAST MATRIX of PRINCIPLE BALANCES        ############################
### This packege is intented to group otus using Compositional Otu Matrix. 
#### Compositional Otu Matrix is - otus vs. samples - a matrix that sum of its column is constant. 
#### It means that samples are composed by otus.     
### This package uses function from minicomposition.R by JJ Egozcue.
### This version was completed in July 2018 by Aslı Boyraz
########################################################################################
# This mini-library contains the following functions:
#### SECTION 1 :Creating New Variation Matrix ####################
# start(X,numberofgroup) : run function on X compositional otu matrix to group otus.
#                          Outputs: Y (Contrast matrix of Groups) and OTUsGroups.
# processColumns(M,numberofgroup) : Process columns of Contrast Matrix of Principle Balances
#                                    to group otus.
# Case_Group_FF(M,colno) : In this situation,row with both + and - 
#                            values have not assigned any group.
# Case_Group_TF(M,colno) : In this situation, row with + values assigned
#                          to a group but rows with - values have not.
#                          assigned any group.
#	ChangeColumn(M,colToGoback) : Changes current processed column. 
#	
#	CheckGroupsPOSorNEG(M,colno) : Checks for group assignment for rows with + and - values. 
#                      Returns flag wth two element: first is for +, second is for - rows.
#	CheckGroups(POSofNEG,M) : Checks for group assignment for rows with + or - values.
#                           Returns flag with one element.(T or F)
#	CombineAssign(M,coln,PS) : Combines row with +,- or both values and assign a groupname.
# NOTCombineJUSTAssign(rows,group): Combines otus and group names in a matrix (OTUsGroups).
# CheckPrevColsPOSorNEG0(M,colno) : Checks each previous column if for POSorNEG and 
#                                   returns flag as long as previous column number.
#
#### SECTION 2: OTHERS ####################
# CreateGroupedOtusDataMatrix(X,numberofgroup) : Creates New Otu Matrix with Groups.
# createContrastMforGroups : creates contrast matrix for Y (output of start() and Checks its validity.. 
# drawPCAbiplot(X) : Draws 2 biplot using 3 PCs of a Matrix. 
# addTaxaInfo(Data,OTUsGroups) : adds taxa information from Data to OTUsGroups which is
#                                output of start(). 
# Table_TaxasOfGroups<-function(OtuTaxa,level) : Create a table with the frequency of taxa per group.
####################################################################################
###### FUNCTIONS IN SECTION 1 ######################################################
###### function start    ######################################################
### Built Contrast Matrix of Principle Balances of X, Create Dendrogram 
### and process columns on M for grouping.
### Returns Contrast Matrix for Grouped otus of x.
#uses mini_compositions_lib.R by JJ Egozcue and library(compositions).
#-----
# X: Compositional Microbiome Data: samples x otus
# numberofgroup: how many group to divide otus.
#-----
# used in : -
# calls: mPBclustvar (mini_compositions_lib.R) , CoDaDendrogram (compositions), processColumns
#        acomp (vegan)
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
start<-function(X,numberofgroup){
  
  # End of the process, this matrix will shows All otus with groups.
  OTUsGroups<<-matrix(, nrow = 0, ncol = 2)
  # Create group names.
  groupnames<<-rev(paste0("G",1:numberofgroup))
  #G10 G9 G8 .. G1
  
  source("mini_compositions_lib.R")
  ## Built Contrast Matrix of Principle Balances
  V<-mPBclustvar(x=X) #samples x otus
  
  ###Represent the CoDa-Dendrogram.
  library(compositions)
  library(vegan)
 # Xdend<-CoDaDendrogram(acomp(X),V=V,type="lines",
  #                      lwd.tree = 1,lwd.leaf = 1,yaxt="n")
  
  #Subset a new table choosing (numberofgroup-1). Create a Groups column. 
  newV<<-as.data.frame(V[,1:(numberofgroup-1)])
  newV$Group<-NA
  
  #At the end of the process, Y matrix will be our new Contrast matrix for OTU Groups. 
  Y<<-newV
  
  ## Process columns of
  newV<-processColumns(newV,numberofgroup)
  print("*** Y matrix ( new Contrast Matrix for grouping ) is created! ")
  print("*** OTUGroups ( Otus and their group info ) is created! ")
 # return(newV)
}

###### function processColumns ######################################################
### Process columns on M for grouping.
### Returns modified M matrix.
#-----
# M: Contrast Matrix of Principle Balances with
# numberofgroup: how many group to divide otus.
#-----
# used in : start
# calls: CombineAssign , CheckGroupsPOSorNEG, Case_Group_FF, Case_Group_TF
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
processColumns<-function(M,numberofgroup){
  
  for(i in rev(c(1:(numberofgroup-1)))){
    
    if(i==(numberofgroup-1)){
      M<-CombineAssign(M,i,"ALL")
    }else{
      flag<-CheckGroupsPOSorNEG(M,i) ## results T-F or F-F
      if(length(unique(flag))==1 && unique(flag)==FALSE){
        M<-Case_Group_FF(M,i)
      }
      if(length(unique(flag))>1){
        M<-Case_Group_TF(M,i)
      }
    }
  }
  return(M)
}

###### function Case_Group_FF ######################################################
### Fires when : flag<-CheckGroupsPOSorNEG(M,colno) is False-False
### In this situation,  row with both + and - values have not assigned any group.
### This function checks previous columns of current column whether any matching
### values. If yes, then that previous column is processed. If no, then current
### column is passed for furthter steps.
### Returns modified M matrix.
#-----
# M: Contrast Matrix of Principle Balances
# colno: the processed column no
#-----
# used in : processColumns
# calls: CheckGroupsPOSorNEG , CheckPrevColsPOSorNEG0, ChangeColumn, 
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
Case_Group_FF<-function(M,colno){
   
  #False - False : check prevs and if no false(matching column) PASS this column! 
  flag<-CheckGroupsPOSorNEG(M,colno) # False False
  
  # if(P == FALSE && N == FALSE)
   if(length(unique(flag))==1){
   #Check prevs. cols
    flags<-CheckPrevColsPOSorNEG0(M,colno)
    print(flags)
	# Look for "False" if exist! Because it points the previous column with matching valaue.
    if(!is.na(which(flags[,1]==FALSE)[1])){ 
      # find matching col
      matchingcol<-which(flags[,1]==FALSE)[1] + colno 
	  cat("A matching column found: Column:",matchingcol) 
      M<-ChangeColumn(M,matchingcol)
      
    }else{
       print("PASS this COlumn for further steps!") 
    }
    
  }
  return(M)
}


###### function Case_Group_TF ####################################
### Fires when : flag<-CheckGroupsPOSorNEG(M,colno) is True-False
### In this situation, row with + values assigned to a group 
### but rows with - values have not assigned any group.
### This function checks previous columns of current column whether 
### any matching values and decide to combine rows and assign group 
### or pass the column for furthter steps.
### Returns modified M matrix.
#-----
# M: Contrast Matrix of Principle Balances
# colno: the processed column no
#-----
# used in : processColumns
# calls: CheckGroupsPOSorNEG , ChangeColumn, 
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
Case_Group_TF<-function(M,colno){
  #TRUE FALSE / FALSE TRUE: for true nothing to do. Check previous cols for False one.
  flag<-CheckGroupsPOSorNEG(M,colno) # True False
	
  POS<-which(M[,colno]>0)
  NEG<-which(M[,colno]<0)
  
  ff<-CheckPrevColsPOSorNEG0(M,colno)

  if(!flag[2]){
    cat("-> Processing NEG..")
    
    f<-ff[,2]
    print(f)
    
    if ( length(unique(f)) == 1 ){ # if all values are true (means all prev cols are 0s.)
      cat("Prev cols are all 0. So just Combine!")
      M<-CombineAssign(M,colno,"NEG") 
     }else { ## eger önceki kolonlarda çakışan var demektir.
        # find matching col
       print("There is a match for prev col (!)")
        matchingcol<-which(f==FALSE) + colno 
        M<-ChangeColumn(M,matchingcol)
     }
    
	 ## THIS PART NOT IN USE ACTUALLY BUT IN ANY CASE IT IS HERE..
     }else if(!flag[1]){
	  cat("-> Processing POS..")
       
       f<-ff[,1]
       print(f)
       if ( uniques(f) == TRUE ){
      cat("Combine POS!")
      #  POSorNEG<-deparse(substitute(POS)) # convert to String.
      M<-CombineAssign(M,colno,"POS")
      }else{
        # find matching col
        matchingcol<-which(f==FALSE) + colno 
        M<-ChangeColumn(M,matchingcol)
      } 
  }
  
  return(M)

}


###### function ChangeColumn ####################################
### calls CombineAssign for the argument colToGoback and M.
### returns modified M matrix.
#-----
# M: Contrast Matrix of Principle Balances
# colToGoback: the processed column no
#-----
# used in : Case_Group_FF and Case_Group_TF
# calls: CombineAssign
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
ChangeColumn<-function(M,colToGoback){
  print(paste("column changed!:",colToGoback))
  temp<-CombineAssign(M,colToGoback,"ALL")
 return(temp)
}

###### function CheckPrevColsPOSorNEG0 ####################################
### Check if previous columns for the rows with positive values and 
### negative values are all 0's.  
### Returns a flag with two columns: 
### -> First column is flags for rows with positive values. 
### -> Second column is flags for rows with negative values.
### Returns a TRUE column if previous columns only contain 0's.  
#-----
# M: Contrast Matrix of Principle Balances
# colno: the processed column no
#-----
# used in : Case_Group_FF
###################################################################
# ASLI BOYRAZ, July 2018
######################################################

CheckPrevColsPOSorNEG0<-function(M,colno){
pos_prevflag<-NULL
neg_prevflag<-NULL

POS<-which(M[,colno]>0)
NEG<-which(M[,colno]<0)

collist<-c((colno+1):(ncol(M)-1))

for(i in collist) {
  print(paste0("Checking for Col.",i))
  if (length(unique(M[POS,i])) < 2 )  { # POS rows on colno is All 0's
    pos_prevflag<- c(pos_prevflag,TRUE) #If all prev cols 0.
  }else{
    pos_prevflag<- c(pos_prevflag,FALSE) ##need to change column
  }
  
}
for(i in collist) {

  if (length(unique(M[NEG,i])) < 2 )  { # NEG rows on colno is All 0's
    neg_prevflag<- c(neg_prevflag,TRUE) #If all prev cols 0.
  }else{
    neg_prevflag<- c(neg_prevflag,FALSE) ##need to change column
  }
}
prevflag<-cbind(pos_prevflag,neg_prevflag)
return(prevflag)
}




###### function CheckGroupsPOSorNEG ####################################
### Check if rows with positive values and negative values has assigned 
### any group previously.
### Returns a flag with two values: First is for to positive and second 
### is for negative values. 
#-----
# M: Contrast Matrix of Principle Balances
# colno: the processed column no
#-----
# used in : Case_Group_TF
# calls : CheckGroups
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
CheckGroupsPOSorNEG<-function(M,colno){
  flag<-NULL
  POS<-which(M[,colno]>0)
  NEG<-which(M[,colno]<0)
  flag<-c(flag,CheckGroups(POS,M)) # T
  flag<-c(flag,CheckGroups(NEG,M)) # F 
  return(flag)
}


###### function CombineAssign ###################################
### According to PN argument, combines values and assigns a group 
### name for the row.
### Returns modified M matrix.  
#-----
# M: Contrast Matrix of Principle Balances
# colno: the processed column no
# PN: "POS" , "NEG" or "ALL"
#   if POS only for rows with positive values.
#	  if NEG only for rows with negative values.
#	  if ALL both for rows with positive and negative values.
#-----
# used in Case_Group_TF
##############################################################
# ASLI BOYRAZ, July 2018
######################################################
CombineAssign<-function(M,colno,PN){
  
  cat("----->COmbining ",PN)
  print("-----")
  
  if(PN=="POS"){
  POS<-which(M[,colno]>0)
  M<-M[-POS[-1],]
  M[which(M[,colno]>0),"Group"]<-groupnames[1]
  
  NOTCombineJUSTAssign(POS,groupnames[1])
  groupnames<<-groupnames[-1]
  Y<<-M
  
  }else if(PN=="NEG"){
  NEG<-which(M[,colno]<0)
  M<-M[-NEG[-1],]
  M[which(M[,colno]<0) ,"Group"]<-groupnames[1]
  
  NOTCombineJUSTAssign(NEG,groupnames[1])
  groupnames<<-groupnames[-1]
  Y<<-M
  
  }else if (PN=="ALL"){
  POS<-which(M[,colno]>0)
  M<-M[-POS[-1],]
  M[which(M[,colno]>0),"Group"]<-groupnames[1]
  
  NOTCombineJUSTAssign(POS,groupnames[1])
  groupnames<<-groupnames[-1]
  Y<<-M
  
  NEG<-which(M[,colno]<0)
  M<-M[-NEG[-1],]
  M[which(M[,colno]<0) ,"Group"]<-groupnames[1]
  
  NOTCombineJUSTAssign(NEG,groupnames[1])
  groupnames<<-groupnames[-1]
  Y<<-M
  }
  
  return(M)
}


###### function NOTCombineJUSTAssign ###################################
### Combines otu groups and assigned groups in OTUsGroups matrix defined in start(). 
#-----
# used in CombineAssign
##############################################################
# ASLI BOYRAZ, July 2018
######################################################
NOTCombineJUSTAssign<-function(rows,group){
  print("----->NOCOmbiningJUSTAssigning ")
  temp<-cbind(rownames(Y[rows,]),group)
  OTUsGroups<<-rbind(OTUsGroups,temp)
}

###### function CheckGroups ########################################
### Checks rows in a column of M whether is assigned a group or not.
### Returns a flag.
#-----
# M: Contrast Matrix of Principle Balances
# POSofNEG: Row_ids with Positive or Negative values in the column
#-----
# used in Case_Group_FF and CombineAssign
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
CheckGroups<-function(POSorNEG,M){
  flag<- NULL
  if (unique(is.na(M[POSorNEG,]$Group))==TRUE) {   
    flag<-c(flag,FALSE)
    print("No assigned group") 
  }else {
    print("Assigned Group:")
    print(unique(M[POSorNEG,]$Group))
    flag<-c(flag,TRUE)
  }
  return(flag)
}

###############################################################################
###############################################################################
########################## FUNCTIONS in Section 2 ##############################
###### function drawPCAbiplot ####################################
### Draws 2 biplot using 2 PC of New Otu Matrix with Groups. 
#-----
# X: Dimension Reduced Compositional Otu Matrix. (samples x otus)
#-----
# used in : -
# calls: mclr (from mini_composition_lib.R) , prcomp, biplot (from stats)
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
drawPCAbiplot<-function(X){
  library(stats)
  ##Draws two biplot using 3 PC.
  X.pc <- prcomp(mclr(X))
  X.mvar <- sum(X.pc$sdev^2)
  PC1 <- paste("PC1: ", round(sum(X.pc$sdev[1]^2)/X.mvar, 3))
  PC2 <- paste("PC2: ", round(sum(X.pc$sdev[2]^2)/X.mvar, 3))
  PC3 <- paste("PC3: ", round(sum(X.pc$sdev[3]^2)/X.mvar, 3))
  biplot(X.pc, choices = 1:2, scale=1,expand =1,xlab=PC1, ylab=PC2,col=c("gray","black"))
  biplot(X.pc, choices = 2:3,scale=1,expand =1,xlab=PC2, ylab=PC3,col=c("gray","black"))
}

###### function CreateGroupedOtusDataMatrix ########################################
### Creates New Otu Matrix with Groups. 
#-----
# X: Compositional Otu Matrix.
# numberofgroup : number of groups.
#-----
# used in : -
# calls: createContrastMforGroups
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
CreateGroupedOtusDataMatrix<-function(X,numberofgroup){
  
  ## Built Contrast Matrix of Principle Balances: V is found using variation tree.
  V<-mPBclustvar(x=X) #samples x otus
  balances<-milr(X,V) #X:samplexotus - milr() assign ilr coordinates
  minibal<-balances[,1:(numberofgroup-1)]

  Y.contrast<-createContrastMforGroups(Y)
  
  ## X.new is 29 x numberofgroups (reduced dimension!)
  X.new<-milrInv(minibal,Y.contrast) # (mileInv from minicomposition_lib.R)
  return(X.new)
}


###### function createContrastMforGroups ########################################
### Creates Contrast matrix for Y and Checks it. 
### Returns a Y.contrast
#-----
# Y: The basis to define new Groups; Output of start().
#-----
# used in : -
# calls : mbuildcontrast (mini_composition_lib.R) , 
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
createContrastMforGroups<-function(Y){
  
  print("Creating Y.contrast Matrix..")
  
  #Sign matrix of Y will be used as basis to built contrast matrix for groups.
  Y.sign<-sign(Y[-numberofgroup])
  rownames(Y.sign)<-Y$Group
  Y.contrast<-mbuildcontrast(Y.sign)
  
  #Check Y.constrast matrix if calculation is true!
  # t(Y.contrast) * Y.contrast should be Identity matrix.
  print("Checking Y.contrast Matrix..")
  print(all(diag(x = 1, ncol(Y.sign), ncol(Y.sign))== round(t(Y.contrast)%*%Y.contrast,2)))
  # Y.contrast * t(Y.contrast) should be symetric.
  #print(unique(round(rowSums(round(Y.contrast%*%t(Y.contrast),2)),2)) == 0)
  print(isSymmetric(round(Y.contrast%*%t(Y.contrast),2)))
  
  return(Y.contrast)
}

###### function addTaxaInfo ########################################
### Adds Taxa info to OTUsGroups MAtrix that is output of start(). 
### Returns a OTUgrp.taxa
#-----
# Data : Initial Data
# OTUsGroups: MAtrix with otus and assigned groups. It is output of start().
#-----
# used in : -
# calls : -
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
addTaxaInfo<-function(Data,OTUsGroups){
  #Data'nın son kolonu otuid önceki 6 kolon tax info
  colnames(OTUsGroups)<-c("otuid","group")
  otulist<-OTUsGroups[,1]
  a<- length(Data)
  Data.filtered<-Data[otulist,c(a,(a-7):(a-1))]
  Data.filtered<-Data[otulist,c(a,(a-7):(a-1))] #otu id and taxa info columns
  OTUgrp.taxa<-merge(OTUsGroups,Data.filtered, by="otuid")
  return(OTUgrp.taxa)
}
###### function addTaxaInfo ########################################
### Create a table with the frequency of taxa per group. Taxa info to OTUsGroups MAtrix that is output of start(). 
### Returns a OTUgrp.taxa
#-----
# level : "Kingdom", "Phylum", "Class" ,"Order" , "Family" , "Genus", "Species"
# OTUsGroups: Matrix with otus,groups and taxonomic info. It is output of addTaxaInfo().
#-----
# used in : -
# calls : - 
###################################################################
# ASLI BOYRAZ, July 2018
######################################################
Table_TaxasOfGroups<-function(OtuTaxa,level){
  library(reshape2)
  temp<-OTUgrp.taxa[,c("group",level)]
  temp<-as.data.frame(table(temp))
  temp$level <- rep(1:10,times = 13)
  # formula is : level ~ group
  tbl<-dcast(data = temp,formula = as.formula(paste(level, "~ group ")), fun.aggregate = sum, value.var = "Freq")
  order<-paste0("G",1:10)
  tbl<-tbl[,c(level,order)]
  return(tbl)
}

