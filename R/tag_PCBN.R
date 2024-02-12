
# Creates a simple string node|cond_set
paste_margin <- function(node, cond_set)
{
  if (length(cond_set) == 0){
    return (node)
  }
  tag = paste0("l_", node, " | ", paste(cond_set, collapse=','), "_r")

  return(tag)
}

# Turns string node|cond_set back to its vector notation
unpaste_margin <- function(tag){
  v = strsplit(tag,split="|", fixed = TRUE)[[1]][1]
  cond_set = strsplit(tag,split="|", fixed = TRUE)[[1]][2]
  cond_set = strsplit(cond_set,split=",", fixed = TRUE)[[1]]
  return(list(v = v, cond_set=cond_set))
}

