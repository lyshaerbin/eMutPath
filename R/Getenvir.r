#'@title Get the environment datasets
#'@export
#'@aliases
initialize_data<-function(){
  utils::data("envData",package="EdgeticDriver")
}

Getenvir<-function(envData){

  if(!exists("envData")) initialize_data()
  return(get(envData,envir="envData"))

}
