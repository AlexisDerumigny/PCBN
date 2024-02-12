
# Creates a simple string node|cond_set
paste_margin <- function(node, cond_set)
{
  if (length(cond_set) == 0){
    return (node)
  }
  tag = paste0(node, " | ", paste(cond_set, collapse=', '))

  return(tag)
}


# Turns string node|cond_set back to its vector notation
unpaste_margin <- function(tag){
  v = strsplit(tag,split="|", fixed = TRUE)[[1]][1]
  cond_set = strsplit(tag,split="|", fixed = TRUE)[[1]][2]
  cond_set = strsplit(cond_set,split=",", fixed = TRUE)[[1]]
  return(list(v = v, cond_set=cond_set))
}


print_key_keychain <- function(key){
  if ( identical(names(key), c("margins", "cond")) ){
    if (length(key[["margins"]]) != 2){
      stop("There should be exactly 2 margins, whereas there are ",
           length(key[["margins"]]), "here.")
    }
    text = paste_margin(node = paste0(key[["margins"]][1], ", ",
                                      key[["margins"]][2]),
                        cond_set = key[["cond"]] )
    return ( text )
  } else if ( identical (names(key), c("margin", "cond")) ){
    if (length(key[["margin"]]) != 1){
      stop("There should be exactly 1 margin, whereas there are ",
           length(key[["margin"]]), "here.")
    }
    text = paste_margin(node = key[["margin"]][1],
                        cond_set = key[["cond"]] )
    return ( text )
  } else {
    stop(" Unknown key, with names: ", names(key) )
  }
}

