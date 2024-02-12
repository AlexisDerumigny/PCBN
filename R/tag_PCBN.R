
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


make_and_store_keyCopula <- function(v, w, cond, v_key, w_key, e)
{
  key_keychain = list(margins = c(v, w), cond = cond)

  # The copula key is just the (ordered) list of the two keys of the
  # (conditional) margins, with its name
  copula_key = list(name = print_key_keychain(key_keychain),
                    margin1 = v_key, margin2 = w_key)

  e$keychain[[key_keychain]] = copula_key
  return (copula_key)
}

make_and_store_keyMargin <- function(v, cond, copula_key, e)
{
  key_keychain = list(margin = v, cond = cond)

  # The margin key is the list of its name, the name of the
  # (conditional) margin 'v' and the copula_key.
  margin_key = list(name = print_key_keychain(key_keychain),
                    margin = v, cond = cond, copula = copula_key)

  e$keychain[[key_keychain]] = margin_key
  return (margin_key)
}


