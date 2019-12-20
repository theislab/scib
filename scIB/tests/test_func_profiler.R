source("~/helmholtz_munich/benchmarking_data_integration/Benchmarking_data_integration/scIB/R_funcs.R")

# load parallel package
require(parallel)
 
# define function to test whether an number is prime
is_prime <- function(num)
{
    # if input equals 2 or 3, then we know it's prime
    if(num == 2 | num == 3) 
      return(TRUE)
    # if input equals 1, then we know it's not prime
    if(num == 1) 
      return(FALSE)
   
    # else if num is greater than 2
    # and divisible by 2, then can't be even
    if(num %% 2 == 0) 
      return(FALSE)
   
    # else use algorithm to figure out
    # what factors, if any, input has
     
    # get square root of num, rounded down
    root <- floor(sqrt(num))
     
    # try to divide each odd number up to root
    # into num; if any leave a remainder of zero,
    # then we know num is not prime
    for(elt in seq(5,root))
    {
        if (num %% elt == 0)
          return(FALSE)
       
    }
    # otherwise, num has no divisors except 1 and itself
    # thus, num must be prime
    return(TRUE)
   
}

test_prime = function(samples, n_threads=3){
	   # create cluster object
	   cl <- makeCluster(n_threads)

	   # test each number in sample_numbers for primality
	   results <- parSapply(cl , sample_numbers , is_prime)
	   
	   # close
	   stopCluster(cl)
}

# get random sample of 10 million integers from integers between 
# 1 and 100 million
# set seed so the random sample will be the same every time
set.seed(2)
sample_numbers <- sample(100000000, 10000000)


print("starting function call 1")
start = proc.time()
res = func_profiler(test_prime(sample_numbers, 1))
end = proc.time()
print("test with 1 core")
print(paste("time:", res$time))
print(paste("memory use:", res$mem))
print(paste("proc time:", (end - start)))

print("")
Sys.sleep(1)
print("")

print("starting function call 2")
start = proc.time()
res = func_profiler(test_prime(sample_numbers, 2))
end = proc.time()
print("test with 2 cores")
print(paste("time:", res$time))
print(paste("memory use:", res$mem))
print(paste("proc time:", (end - start)))

print("")
Sys.sleep(1)
print("")

print("starting function call 3")
start = proc.time()
res = func_profiler(test_prime(sample_numbers, 3))
end = proc.time()
print("test with 3 cores")
print(paste("time:", res$time))
print(paste("memory use:", res$mem))
print(paste("proc time:", (end - start)))


# Results:
# 1 proc:
# - time: 7.24
# - memory: 1730.2 MB
#
# 2 procs:
# - time: 6.8
# - memory: 2106.3 MB
#
# 3 procs: 
# - time: 5.94
# - memory: 2076.5 MB
#
# Note 1: this seems to also measure the data passed to the function.
# Note 2: Real runtime seemed a lot longer than the profiled 7s.
