library(deSolve)
# Simulate netowrk dynamics using the SQUAD method for standardized qualitative analysis of dynamical networks
# Analyze the effect of mutating interactions and nodes on the network

setwd(".")

mutateNodes<-function(node.names, value, rules){
  rules.tmp<-rules
  mutant.rules<-lapply(seq_along(node.names), function(node){
        rules.tmp[node]<-value
        gsub(paste0("\\b",as.character(node.names[node]),"\\b"),value,rules.tmp)
    })
  names(mutant.rules)<-node.names
  return(mutant.rules)
}

parseRules<-function(rules){
  sapply(rules, function(rule){
  parse(text=rule)})
}



# evalRules<-function(expressions){
#               sapply(expressions, function(x){
#               eval(x)})}

squad<-function(x,w,h,g){
((-exp(0.5*h) + exp(-h*(w-0.5))) / ((1-exp(0.5*h)) * (1+exp(-h*(w-0.5))))) - (g*x)
}


network<-function(t, state, parameters) {
  
  with(as.list(c(state, parameters)),{
  
    eval.pars<-sapply(parameters, function(x){eval(x)}) 
    derivatives<-lapply(seq_along(state), function(x){
		  unname(squad(state[x], eval.pars[x],h=50, g=1.0))
      })

    #Model equations using SQUAD
    
    #Results are returned as a list
    return(list(derivatives))})}

simulateMutant<-function(mutant.rules, variables,deriv.vars,node.value, i, x,time, parameters, filename, current.node){
mutant.rules<-unlist(mutant.rules)

  cat("Node:", current.node,"  Realization", i, "\n")  
  #Assing a random number to each variable
  lapply(names(variables), function(x){ 
	if(x==current.node){assign(x, value=node.value, .GlobalEnv)}
	else{assign(x,value=runif(1,0,1), .GlobalEnv)}
	})
	
  #initial state
  state<-c(mget(names(variables), envir= .GlobalEnv), recursive=TRUE)
  #Evaluate parsed rules and assign values to general parameters of the ODEs system
  
  #parameters<-c(parameters,mutant.rules)
  parameters<-mutant.rules
  #ODE result (vector with time and node variables)
  out <- ode(y=state,times=time,func=network,parms=parameters, atol=10e-6, rtol=10e-6)
  #Keep only the values of node variables at the final state 
  
  last_row<-out[nrow(out),2:ncol(out)]
  #Save by row and wirte results to a file 
  file<-paste(filename,"_",current.node,".txt")
  
  if(!file.exists(file)){
    write.table(t(last_row),file, sep=",", append=FALSE, quote = FALSE, row.names = FALSE, col.names = names(variables))
  }
  
  else{
    write.table(t(last_row),file, sep=",", append=TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }


}

mutantAnalysis<-function(rules, node.names, network.name, type=c("overexp", "null"), realizations, time,  parameters){
#cat(type)
if(missing(node.names)){
  #node.names<-c("Akt","B9","Bach2","Bcl2","Bcl6","Bcl11b","BCR","Blimp1","CD4","CD8","CD19","Dll1","Ebf1","EOMES","ERK","Flt3","Flt3L","Foxp3","FR4","GATA3","Gfi1","Gzmb","HEB","Helios","Hoxa9","IFNb","IFNbR","IFNg","IFNgR","Ikaros","IL2","IL2R","IL4","IL4R","IL6","IL6R","IL7","IL7R","IL10","IL10R","IL12","IL12R","IL17","IL18","IL18R","IL21","IL21R","IL23","IL23R","IL27","IL27R","IRAK","Irf4","JAK1","JAK3","NFkB","NFAT","Notch1","Pax5","Prf1","PU1","RORgt","Runx1","Runx3","SMAD2","SMAD3","SOCS1","STAT1","STAT3","STAT4","STAT5","STAT6","Tbet","TCF1","TCR","TGFb","TGFbR","ThPOK","TNFa","TNFR2","XBP1")
  stop("Please supply a vector of node names")
}

if(match.arg(type)=="overexp"){
  node.value<-1.0
  filename<-paste(network.name,"Mutant_Analysis_Overexp", sep="_")    
    mutants.list<-mutateNodes(node.names,1, rules = rules)
		
}
if(match.arg(type)=="null"){
  node.value<-0
  filename<-paste(network.name,"Mutant_Analysis_NullExp", sep="_")  
  mutants.list<-mutateNodes(node.names,0, rules = rules)
  
}



if(missing(time)){
  time <- seq(0,100)
}

if(missing(parameters)){
  parameters<-c(h=50,g=1.0)
}

#network.rules<-rules
number.of.nodes<-length(rules)
node.names<-node.names

network.name<-network.name

realizations<-as.integer(realizations)
#Create variables in .Globalenv for each node and its derivative.
deriv.vars<-vector("list", length(node.names))
names(deriv.vars)<-sapply(node.names,function(x){paste0("d",x)})
list2env(deriv.vars, env=.GlobalEnv)


variables<-vector("list", length(node.names))
names(variables)<-node.names
list2env(variables, env=.GlobalEnv)

rules.parsed<-parseRules(rules)


mutants.parsed.rules<-lapply(mutants.list, function(x){parseRules(x)})
names(mutants.parsed.rules)<-node.names
for(i in 1:realizations){
sapply(seq_along(mutants.parsed.rules), function(x){ 
  current.node<-names(mutants.parsed.rules)[x]
  simulateMutant(mutants.parsed.rules[[x]], variables, deriv.vars, node.value, i,x, time, parameters, filename, current.node)})
}
} 
#Usage example with the Lymphopoiesis network from Mendoza and Mendez,2015. The regulatory rules are passed to the function in 
#the form of a vector and the type of perturbation is chosen from "overexp" or "null" to simulate single gain- and loss-of-function mutations, respectively.
#The user has to supply the name of the file to save the simulation results and the number of simulations (realizations) to perform for each mutant. 

#rules<-c("TNFR2","0","max(Pax5,Bcl6)","max(FR4,STAT5)","min(IL21R,Ebf1,1-Blimp1,1-Irf4)","min(Notch1,TCF1)","0","min(Irf4,1-Bach2,1-Bcl6,1-Pax5)","min(max(CD4,Notch1,ThPOK),1-Ebf1,1-Runx3)","min(max(CD8,Notch1,Runx3),1-Ebf1,1-TCR,1-ThPOK)","Pax5","0","min(Runx1,1-Runx3)","min(max(IL27R,Tbet),Runx3)","max(BCR,Flt3)","min(max(Flt3L,Hoxa9),1-Pax5)","1","min(max(SMAD2,SMAD3,Foxp3),1-GATA3,1-RORgt,1-STAT3,1-Tbet)","min(B9,Foxp3)","min(max(GATA3,STAT6,TCF1),1-Foxp3,1-HEB,1-RORgt,1-Runx1,1-Tbet)","Ikaros","EOMES","Notch1","NFkB","0","0","IFNb","min(max(EOMES,IRAK,STAT4,Tbet),1-STAT3)","IFNg","PU1","0","IL2","min(GATA3,1-STAT1)","min(IL4,1-SOCS1)","RORgt","IL6","0","min(IL7,1-TCR)","GATA3","IL10","0","min(IL12,1-STAT6)","RORgt","0","min(IL18,1-STAT6)","0","IL21","0","IL23","0","IL27","IL18R","min(max(Blimp1,NFkB),Ebf1,1-Flt3)","IFNgR","min(max(IL6R,IL7R),1-Notch1)","BCR","TCR","Dll1","max(min(Ebf1,1-Blimp1,1-Flt3,1-Irf4,1-ThPOK),min(Ebf1,STAT5,1-Blimp1,1-Irf4,1-ThPOK))","EOMES","min(1-Gfi1,1-Ikaros)","max(min(SMAD2,STAT3),min(RORgt,1-Foxp3,1-GATA3,1-Tbet))","min(Runx1,1-Notch1,1-Runx3,1-TCR,1-ThPOK)","max(min(CD8,1-CD4,1-STAT5,1-ThPOK),min(CD8,STAT5))","TGFbR","min(TGFbR,1-Akt)","max(STAT1,Tbet)","max(IFNbR,IL27R,JAK1)","max(IL10R,IL21R,IL23R,JAK3)","min(IL12R,1-GATA3)","min(max(IL2R,JAK3),1-SOCS1)","IL4R","min(max(STAT1,Tbet),1-Foxp3,1-GATA3,1-RORgt)","Notch1","0","Foxp3","TGFb","max(min(CD4,1-CD8),Foxp3,GATA3,RORgt,TCR,min(Tbet,1-Runx3),ThPOK)","RORgt","TNFa","Blimp1")
#mutantAnalysis(rules, network.name="Net_81Nodes_1000_Realizations", type="overexp", realizations=1000)
