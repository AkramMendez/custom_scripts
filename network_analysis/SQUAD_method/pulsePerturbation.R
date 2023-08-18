#This function allows to simulate multiple perturbations by creating events 
#that can then be applied to the system of equations. 
# 	
#	Usage: The user can provide a vector with the names of nodes to perturb in the form: c(Node_1, Node_2,..., Node_n), 
#	the 'time' or times at which the perturbations will be applied c(time_1, time_2,..., time_n),
#	the 'duration' of the time interval for each perturbation c(duration_1, duration_2,...,duration_n),
# 	the 'type' of perturbation to perform: "on" to simulate the activation of a node (value = 1.0), 
#	"off" to simulate the inactivation of a node (value = 0) 
# 	for a given interval of time and "replace" to simulate any other value.
#	If 'type'=="other", the 'intensity' parameter values for each perturbation must be provided in the range [0,1].

createEvent<-function(nodes, time, duration, type=c("on", "off", "other"), intensity, simulation.time){
	if(missing(nodes))
	stop("Supply the name of a single node or a list of nodes to perturbate")
	if(missing(time))
	 stop("Supply a value or vector of values of the time at which the perturbation will be applied")
	 if(missing(duration))
	 stop("Supply the duration or multiple durations by which the perturbations will be applied")
	 if(match.arg(type)=="other" & missing(intensity))
	 stop("Please supply a vector with values in the range [0,1] for each perturbation")
	if(missing(type))
		stop("Supply the type of perturbation to perform")
		
	time.step.size<-(max(simulation.time)-min(simulation.time))/(length(simulation.time)-1)
	adjust.time<-(time/time.step.size)+1
	adjust.duration<-duration/time.step.size
	from.time<-adjust.time
	to.time<-adjust.time+adjust.duration
	#interval<-(duration/time.step.size)+1

	variables<-pulses<-intensities<-NULL


	for(i in 1:length(from.time)){
		pulse<-seq(from.time[i], to.time[i], by=time.step.size)
		pulses<-c(pulses,pulse)	
		variables<-c(variables,rep(nodes[i],length(pulse)))
			if(match.arg(type)=="other" & !missing(intensity)){
			intensities<-c(intensities, rep(intensity[i], length(pulse)))
			}
	}
	
	if(match.arg(type)=="on"){
		intensities<-rep(1.0, length(variables))
	}
	if(match.arg(type)=="off"){
		intensities<-rep(0, length(variables))
	}

	types<-rep("rep", length(variables))

	#perturbation.events<-data.frame(var=variables, time=pulses, value=intensities, method=types, stringsAsFactors=FALSE)
	perturbation.events<-data.frame(var=variables, time=pulses, value=intensities, method=types)

	return(perturbation.events)
}

createEvent(names(initial.state), time = c(5,10,15,20,25), duration = c(3,4,5,6,3), intensity = c(1,1,1,1,1), type = "on")