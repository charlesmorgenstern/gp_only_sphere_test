#####################################################
#####################################################
using Plots
using Printf
#####################################################
#####################################################
function getpaths(n,r)
#get 6 paths with n points starting at radius r
#paths are x,y,z axis

t=collect(range(r,1,n)) #parameter for gps
p=Matrix{Float64}(undef,n,18) #array to hold gps
p[:,:].=0.0 #set to zero

p[:,1]=t #positive x axis
p[:,4]=-t #negative x axis
p[:,8]=t #postive y axis
p[:,11]=-t #negative y axis
p[:,15]=t #positive z axis
p[:,18]=-t #negative z axis

return p
end
#####################################################
#####################################################
function int_vol(n,r)
#approximate volume of the unit sphere with gp only integration
#use 6 paths from getpaths with n points for integration
#start paths on interior sphere of radius r

p=getpaths(n,r) #generate paths

exact=(4/3)*pi #exact volume of unit sphere
intv=(4/3)*pi*r^3 #volume of interior sphere

curv(x,y,z)=1/sqrt(x^2+y^2+z^2) #mean curvature of isosurface at (x,y,z)

h=(1-r)/(n-1) #spacing in gp points

element=(4/6)*pi*r^2 #area of surface element on interior sphere

da=Matrix{Float64}(undef,n,6) #intialize array for da
da[1,:].=element #first entry is surface area of element

for i=1:6 #loop over gps
for j=2:n #loop over points for gp
    int=0
    for k=2:j
    int+=2*curv(p[k,3*i-2],p[k,3*i-1],p[k,3*i])*h #integrate curvature rh rule
    end
    da[j,i]=da[1,i]*exp(int) #take exponential and multiply by element area
end
end

#trapezoidal rule to integrate dA(s) down each path
vol=sum(da[1,:])+sum(da[end,:]) #first and last terms

for i=2:n-1
vol+=sum(da[i,:])*2 #middle terms
end

vol*=h/2 #weight
vol+=intv #add volume of interior sphere where gps start
err=abs(vol-exact)/exact #calculate relative error
return err
end
#######################################################
#######################################################
#####################################################
function errortable()
@printf "\n Approximate area of unit sphere using gp-only integration"
@printf "\n Interior sphere radius where gps start: .05"
@printf "\n Number of gradient paths: 6"
@printf "\n n is number of points on each path"
@printf "\n Error is relative error"
@printf "\n EOC is estimated order of convergence"
@printf "\n ----------------------------------------------------------------"
@printf "\n n           rel. error         EOC"

n=(10,20,40,80,160,320,640,1280,2560,2560*2,2560*4)
er=Array{Float64}(undef,11)
eoc=Array{Float64}(undef,11)

er[1]=int_vol(n[1],.05)
@printf "\n %g           %g          n/a" n[1] er[1]

for i=2:11
er[i]=int_vol(n[i],.05)
eoc[i]=log(er[i-1]/er[i])/log(2)
@printf "\n %g           %g          %g" n[i] er[i] eoc[i]
end

end

