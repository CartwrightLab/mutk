
mom_low <- c(1, 0, 0)
dad_low <- c(1, 1e-4, 0)
cld_low <- c(0, 1, 1e-5)
mom_upp <- c(1-1e-4-1e-8, 1e-4, 1e-8)
dad_upp <- c(1-1e-4-1e-8, 1e-4, 1e-8)

# P(i,j) = P(i|j)
# j -> i
make_mat1 <- function(t) {
    matrix(c(1-t,t,t,1-t),2,2)
}
# a/b -> x/y
make_mat2 <- function(t) {
    m <- make_mat1(t)
    k <- kronecker(m,m)[,c(1,2,4)]
    rbind(k[1,],k[2,]+k[3,],k[4,])
}
# a/b -> x
make_mat3 <- function(t) {
    m <- make_mat1(t)
    cbind(m %*% c(1,0),
          m %*% c(0.5,0.5),
          m %*% c(0,1))
}
# a/b x c/d -> x/y
make_mat4 <- function(t1,t2) {
    m1 <- make_mat3(t1)
    m2 <- make_mat3(t2)
    k <- kronecker(m1,m2)
    rbind(k[1,],k[2,]+k[3,],k[4,])
}

# peel
cld <- t(make_mat4(0,0.01)) %*% cld_low
mom <- t(make_mat2(0.001)) %*% mom_low
dad <- t(make_mat2(0)) %*% dad_low
par <- kronecker(dad*dad_upp,mom*mom_upp)
sum(cld * par)


sum((t(make_mat2(0.001)) %*% c(1,0.1,0.001)) * dad_upp)

