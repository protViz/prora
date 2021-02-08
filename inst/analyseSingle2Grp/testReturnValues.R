args = commandArgs(trailingOnly = TRUE)
what <- args[1]
if (what == "q") {
  print("Quitting")
  q(status = 12)
} else if (what == "s") {
  stop("I am stopping")
}else if (what == "run") {
  print("hello")
}
