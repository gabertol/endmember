#' main function "RECA"
#' q= number of end-members (eg. 3)
#' c5= weighting exponent (suggested = 1)
#' i= no. of iterations (eg 1000)
#' c1 = convexity threshold (e.g. -6)
#' @export


RECA=function (x)
{
  library(compositions)
  library(e1071)
  library(nnls)

  #_______________________________________________________________________

  x <- as.matrix(x)
  x <- t(apply(x, 1, norm1))		#__normalize x to a constant row sum
  x_max <- apply(x, 2, max)

  p <- ncol(x)

  EM_dim=10
  if(ncol(x)<EM_dim)EM_dim=ncol(x)

  #________________estimation number of End-Member (EM)______________

  r_fit <- matrix(nrow = p, ncol = EM_dim)
  w <- x
  dec=svd(w)
  B_ini <- dec$v
  ev <- dec$d
  x.pr=list()
  for (l in seq(2, EM_dim, 1)) {
    B_red <- B_ini[, 1:l]
    m_ini <- w %*% B_red
    xprime <- m_ini %*% t(B_red)
    r_fit[, l] <- r2(x, xprime)$R2
    x.pr[[l]]<-xprime
  }
  plot.new
  par(mfrow = c(2, 1))
  numbers <- seq(2, EM_dim, 1)
  numbers <- as.character(numbers)
  matplot(r_fit[, 2:EM_dim], type = "o", lty = 2, cex = 0.8, pch = numbers,
          col = rainbow(length(numbers)), ylab = "R2")
  matplot(cumsum(ev)/sum(ev), main = "normalized eigenvalues",
          type = "o", pch = 1, cex = 0.6)

  #________________________start model (clustering)___________________________________________


  q <- as.numeric(readline(paste("number of endmembers:  ")))							#__manual Input number of EM
  c1 <- as.numeric(readline(paste("please enter convexity threshold (e.g. -6):  ")))	#__manual input of convexity threshold


  repeat {								#__repeat-loop for verification of the start model

    B_ini <- cmeans(x.pr[[q]], q, iter.max = 100, verbose = TRUE,
                    dist = "euclidean", method = "cmeans", m = 1.1, rate.par = NULL)$centers

    #B_ini <- kmeans(x,centers=q,iter.max=100,nstart=5)$centers

    plot.new
    par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = FALSE, par(mfrow = c(1,1)))
    matplot(t(B_ini), type = "l", lty = 1, lwd = 1.5)
    print(paste("rank of initial B: ", qr(B_ini)$rank))
    a <- as.numeric(readline(paste("is the start model plausible? yes = 1; no = 0 :  ")))
    if (a == 1) break
  }


  B <- B_ini
  m_ini <- t(solve(B %*% t(B)) %*% B %*% t(x.pr[[q]]))
  m <- m_ini

  con_err <- c_err(m_ini)$con_err

  print(paste("convexity error: ", con_err))

  xprime <- m_ini %*% B

  if (con_err <= c1){             #check for initial CE
    list(x_mod = xprime,
         B_mod = t(B_ini),
         m_mod = m_ini,
         Convexity_error = con_err,
         ev=ev,
         c1=c1)
  }


  if (con_err <= c1)print("the convexity error of the initial model is already below the treshold!!!")

  if (con_err > c1){

    B_list = list()
    m_list = list()
    Convexity_error = list()

    B_list[[1]] <- B
    m_list[[1]] <- m
    Convexity_error[[1]] <- con_err

    i <- as.numeric(readline(paste("number of iterations:     ")))				#input of max. interation number
    c5 <- as.numeric(readline(paste("please enter weighting exponent c5: ")))	#input of weighting exponent

    exp.err=list()
    #e <- 0
    #count <- c(1:i)

    ########################################################

    for (j in c(1:i)) {
      print("__________________________________________________________________") #__ the iteration cycle

      print(paste("iteration no.: ", j))#count[j-1]))
      if (con_err <= c1) break

      T <- tcalc(m_list[[j]], c5)$T					#calculation of T-matrix, see Appendix Weltje 1997
      B_neu <- T %*% B_list[[j]]
      #B_neu[B_neu < 0] <- 0
      #B_neu <- t(apply(B_neu, 1, norm1))

      m_neu <- t(solve(B_neu %*% t(B_neu)) %*% B_neu %*% t(xprime))
      con_err_2 <- c_err(m_neu)$con_err
      print(paste("Convextiy error:  ", con_err_2))
      con_err <- con_err_2


      Convexity_error[[j + 1]] <- con_err_2
      B_list[[j + 1]] <- B_neu
      m_list[[j + 1]] <- m_neu

      if (min(m_neu) >= 0) break
      if (con_err <= c1)   break
    }




    ################################################################################################
    print("----------------------SUMMARY------------------------")


    if(con_err>=c1)print(">>>> the convexity error threshold was not reached <<<<")

    Convexity_error <- do.call(rbind, Convexity_error)


    k <- length(B_list)
    if(con_err>=c1) k=which.min(Convexity_error)

    B_mod <- B_list[[k]]
    B_mod[B_mod < 0] <- 0
    B_mod <- t(apply(B_mod, 1, norm1))
    m_mod <- m_list[[k]]

    print(paste("convexity error threshold:  ",c1))
    print(paste("weighting exponent:  ",c5))
    print(paste("no: of iterations: ", k-1))
    print(paste("minima in M:  ", min(m_mod)))
    print(paste("negative mean in M:  ", mean(m_mod[m_mod < 0])))
    print(paste("negative values in M :", length(m_mod[m_mod <0]), " of:", length(m_mod), "| in % :",
                (length(m_mod[m_mod < 0]) * 100/length(m_mod))))
    print(paste("Rank of B:  ", qr(B_mod)$rank))
    print(paste("final convexity error  ", con_err_2))

    neg.min=min(m_mod)
    neg.mean=mean(m_mod[m_mod < 0])
    neg.perc=length(m_mod[m_mod < 0]) * 100/length(m_mod)

    #removing neg. values from m_mod
    for (k in c(1:nrow(x))) {if (min(m_mod[k, ]) < 0) {m_mod[k, ] <- nnls(t(B_mod), x[k, ])$x
    }
    }

    m_mod <- t(apply(m_mod, 1, norm1))
    x_mod <- m_mod %*% B_mod

    plot.new
    par(mfrow = c(2, 3))
    matplot(t(B_mod), type = "l", lty = 2, col = "darkgrey")
    matplot(t(B_ini), type = "l", lty = 1, ylab = "B_ini", add = TRUE)
    matplot(t(x), type = "l", lty = 1, cex = 0.7, col = "lightgrey",
            ylim = c(0, max(B_mod)))

    matplot(t(B_mod), type = "o", add = TRUE, lty = 1, cex = 0.7)
    abline(h = 0, lty = 2, col = "grey")


    matplot(Convexity_error, type = "l", xlab = "no. of iterations")

    plot(diag(cor(x, xprime)), type = "o", ylim = c(0, 1))
    points(diag(cor(x, x_mod)), type = "o", pch = 4, col = "blue",
           ylim = c(0, 1))

    matplot(m_ini, type = "l", lty = 1)
    abline(h = 0, lty = 2, col = "grey")
    matplot(m_mod, type = "l", lty = 1)
    abline(h = 0, lty = 2, col = "grey")
    par(mfrow=c(1,1))
  }
  list(x_mod = x_mod, B_mod = B_mod,B_ini=B_ini, m_mod = m_mod,m_ini=m_ini, ev = ev,
       CE = Convexity_error,neg.min=neg.min,neg.mean=neg.mean,neg.perc=neg.perc,c1=c1,c5=c5)



}


RECA=function (x)
{
  library(compositions)
  library(e1071)
  library(nnls)

  #_______________________________________________________________________

  x <- as.matrix(x)
  x <- t(apply(x, 1, norm1))		#__normalize x to a constant row sum
  x_max <- apply(x, 2, max)

  p <- ncol(x)

  EM_dim=10
  if(ncol(x)<EM_dim)EM_dim=ncol(x)

  #________________estimation number of End-Member (EM)______________

  r_fit <- matrix(nrow = p, ncol = EM_dim)
  w <- x
  dec=svd(w)
  B_ini <- dec$v
  ev <- dec$d
  x.pr=list()
  for (l in seq(2, EM_dim, 1)) {
    B_red <- B_ini[, 1:l]
    m_ini <- w %*% B_red
    xprime <- m_ini %*% t(B_red)
    r_fit[, l] <- r2(x, xprime)$R2
    x.pr[[l]]<-xprime
  }
  plot.new
  par(mfrow = c(2, 1))
  numbers <- seq(2, EM_dim, 1)
  numbers <- as.character(numbers)
  matplot(r_fit[, 2:EM_dim], type = "o", lty = 2, cex = 0.8, pch = numbers,
          col = rainbow(length(numbers)), ylab = "R2")
  matplot(cumsum(ev)/sum(ev), main = "normalized eigenvalues",
          type = "o", pch = 1, cex = 0.6)

  #________________________start model (clustering)___________________________________________


  q <- as.numeric(readline(paste("number of endmembers:  ")))							#__manual Input number of EM
  c1 <- as.numeric(readline(paste("please enter convexity threshold (e.g. -6):  ")))	#__manual input of convexity threshold


  repeat {								#__repeat-loop for verification of the start model

    B_ini <- cmeans(x.pr[[q]], q, iter.max = 100, verbose = TRUE,
                    dist = "euclidean", method = "cmeans", m = 1.1, rate.par = NULL)$centers

    #B_ini <- kmeans(x,centers=q,iter.max=100,nstart=5)$centers

    plot.new
    par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = FALSE, par(mfrow = c(1,1)))
    matplot(t(B_ini), type = "l", lty = 1, lwd = 1.5)
    print(paste("rank of initial B: ", qr(B_ini)$rank))
    a <- as.numeric(readline(paste("is the start model plausible? yes = 1; no = 0 :  ")))
    if (a == 1) break
  }


  B <- B_ini
  m_ini <- t(solve(B %*% t(B)) %*% B %*% t(x.pr[[q]]))
  m <- m_ini

  con_err <- c_err(m_ini)$con_err

  print(paste("convexity error: ", con_err))

  xprime <- m_ini %*% B

  if (con_err <= c1){             #check for initial CE
    list(x_mod = xprime,
         B_mod = t(B_ini),
         m_mod = m_ini,
         Convexity_error = con_err,
         ev=ev,
         c1=c1)
  }


  if (con_err <= c1)print("the convexity error of the initial model is already below the treshold!!!")

  if (con_err > c1){

    B_list = list()
    m_list = list()
    Convexity_error = list()

    B_list[[1]] <- B
    m_list[[1]] <- m
    Convexity_error[[1]] <- con_err

    i <- as.numeric(readline(paste("number of iterations:     ")))				#input of max. interation number
    c5 <- as.numeric(readline(paste("please enter weighting exponent c5: ")))	#input of weighting exponent

    exp.err=list()
    #e <- 0
    #count <- c(1:i)

    ########################################################

    for (j in c(1:i)) {
      print("__________________________________________________________________") #__ the iteration cycle

      print(paste("iteration no.: ", j))#count[j-1]))
      if (con_err <= c1) break

      T <- tcalc(m_list[[j]], c5)$T					#calculation of T-matrix, see Appendix Weltje 1997
      B_neu <- T %*% B_list[[j]]
      #B_neu[B_neu < 0] <- 0
      #B_neu <- t(apply(B_neu, 1, norm1))

      m_neu <- t(solve(B_neu %*% t(B_neu)) %*% B_neu %*% t(xprime))
      con_err_2 <- c_err(m_neu)$con_err
      print(paste("Convextiy error:  ", con_err_2))
      con_err <- con_err_2


      Convexity_error[[j + 1]] <- con_err_2
      B_list[[j + 1]] <- B_neu
      m_list[[j + 1]] <- m_neu

      if (min(m_neu) >= 0) break
      if (con_err <= c1)   break
    }




    ################################################################################################
    print("----------------------SUMMARY------------------------")


    if(con_err>=c1)print(">>>> the convexity error threshold was not reached <<<<")

    Convexity_error <- do.call(rbind, Convexity_error)


    k <- length(B_list)
    if(con_err>=c1) k=which.min(Convexity_error)

    B_mod <- B_list[[k]]
    B_mod[B_mod < 0] <- 0
    B_mod <- t(apply(B_mod, 1, norm1))
    m_mod <- m_list[[k]]

    print(paste("convexity error threshold:  ",c1))
    print(paste("weighting exponent:  ",c5))
    print(paste("no: of iterations: ", k-1))
    print(paste("minima in M:  ", min(m_mod)))
    print(paste("negative mean in M:  ", mean(m_mod[m_mod < 0])))
    print(paste("negative values in M :", length(m_mod[m_mod <0]), " of:", length(m_mod), "| in % :",
                (length(m_mod[m_mod < 0]) * 100/length(m_mod))))
    print(paste("Rank of B:  ", qr(B_mod)$rank))
    print(paste("final convexity error  ", con_err_2))

    neg.min=min(m_mod)
    neg.mean=mean(m_mod[m_mod < 0])
    neg.perc=length(m_mod[m_mod < 0]) * 100/length(m_mod)

    #removing neg. values from m_mod
    for (k in c(1:nrow(x))) {if (min(m_mod[k, ]) < 0) {m_mod[k, ] <- nnls(t(B_mod), x[k, ])$x
    }
    }

    m_mod <- t(apply(m_mod, 1, norm1))
    x_mod <- m_mod %*% B_mod

    plot.new
    par(mfrow = c(2, 3))
    matplot(t(B_mod), type = "l", lty = 2, col = "darkgrey")
    matplot(t(B_ini), type = "l", lty = 1, ylab = "B_ini", add = TRUE)
    matplot(t(x), type = "l", lty = 1, cex = 0.7, col = "lightgrey",
            ylim = c(0, max(B_mod)))

    matplot(t(B_mod), type = "o", add = TRUE, lty = 1, cex = 0.7)
    abline(h = 0, lty = 2, col = "grey")


    matplot(Convexity_error, type = "l", xlab = "no. of iterations")

    plot(diag(cor(x, xprime)), type = "o", ylim = c(0, 1))
    points(diag(cor(x, x_mod)), type = "o", pch = 4, col = "blue",
           ylim = c(0, 1))

    matplot(m_ini, type = "l", lty = 1)
    abline(h = 0, lty = 2, col = "grey")
    matplot(m_mod, type = "l", lty = 1)
    abline(h = 0, lty = 2, col = "grey")
    par(mfrow=c(1,1))
  }
  list(x_mod = x_mod, B_mod = B_mod,B_ini=B_ini, m_mod = m_mod,m_ini=m_ini, ev = ev,
       CE = Convexity_error,neg.min=neg.min,neg.mean=neg.mean,neg.perc=neg.perc,c1=c1,c5=c5)



}
