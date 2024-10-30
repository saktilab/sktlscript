#!/usr/bin/env ruby

class Spline
    attr_reader :x0, :x1, :a0, :a1, :a2, :a3, :a4, :a5
    def initialize(x0,x1,a0,a1,a2,a3,a4=0.0, a5=0.0)
        @x0 = x0
        @x1 = x1
        @a0 = a0
        @a1 = a1
        @a2 = a2
        @a3 = a3
        @a4 = a4
        @a5 = a5
    end
    def eval(r, deriv = 0)
        dx = r-@x0
        if (deriv == 0)
            return @a0+@a1*dx+@a2*dx*dx+@a3*dx*dx*dx+@a4*dx*dx*dx*dx+@a5*dx*dx*dx*dx*dx
        elsif (deriv == 1)
            return @a1+2*@a2*dx+3*@a3*dx*dx+4*@a4*dx*dx*dx+5*@a5*dx*dx*dx*dx
        elsif (deriv == 2)
            return 2*@a2+3*2*@a3*dx+4*3*@a4*dx*dx+5*4*@a5*dx*dx*dx
        end
    end
end

class Spline4
    attr_reader :x0, :x1, :a0, :a1, :a2, :a3, :a4
    def initialize(x0,x1,a0,a1,a2,a3,a4)
        @x0 = x0
        @x1 = x1
        @a0 = a0
        @a1 = a1
        @a2 = a2
        @a3 = a3
        @a4 = a4
    end
    def eval(r, deriv = 0)
        dx = r-@x0
        if (deriv == 0)
            return @a0+@a1*dx+@a2*dx*dx+@a3*dx*dx*dx+@a4*dx*dx*dx*dx
        elsif (deriv == 1)
            return @a1+2*@a2*dx+3*@a3*dx*dx+4*@a4*dx*dx*dx
        elsif (deriv == 2)
            return 2*@a2+3*2*@a3*dx+4*3*@a4*dx*dx
        end
    end
end

x0s = []
splines = []
a = 0
b = 0
c = 0
deriv = 0
if (ARGV[1])
    deriv = ARGV[1].to_i
end
File.open(ARGV[0], "r") do |file|
    while line = file.gets
        if (line.include?("Spline")) 
           break 
        end
    end
    splineorder = 3
    if (line.include?("Spline4"))
        splineorder = 4
    end
    line = file.gets
    arr = line.split(" ")
    nspline = arr[0].to_i
    line = file.gets
    a, b, c = line.split(" ").map{|x| x.to_f}
    

    (1..nspline).each do |i|
        line = file.gets
        arr = line.split(" ")
        x0 = arr[0].to_f
        x1 = arr[1].to_f
        
        x0s << x0
#        spline = lambda { |x| (arr[2].to_f)+(arr[3].to_f()*(x-x0))+(arr[4].to_f()*(x-x0)**2)+(arr[5].to_f()*(x-x0)**3)+(arr[6].to_f()*(x-x0)**4) }
#        if ( i == nspline )
#            spline = lambda { |x| (arr[2].to_f)+(arr[3].to_f()*(x-x0))+(arr[4].to_f()*(x-x0)**2)+(arr[5].to_f()*(x-x0)**3)+(arr[6].to_f()*(x-x0)**4)+arr[7].to_f()*(x-x0)**5 }
#        end
        if (splineorder == 3)
            spline = Spline.new(x0, x1, arr[2].to_f, arr[3].to_f, arr[4].to_f, arr[5].to_f, arr[6].to_f, arr[7].to_f)
        elsif (splineorder == 4)
            spline = Spline4.new(x0, x1, arr[2].to_f, arr[3].to_f, arr[4].to_f, arr[5].to_f, arr[6].to_f)
        else
            puts "ERROR"
            exit 2
        end
        splines << spline

#        int = ((x1-x0)/10.0)
#        (0..9).each do |j|
#            puts "%.6f %.6f" % [x0+j*int, spline.call(x0+j*int)]
#        end
        if (i==nspline)
#            puts "%.6f %.6f" % [x1, spline.call(x1)]
            x0s << x1
        end
    end
end 


#a = -2.0*splines[0].a2/splines[0].a1
#c = splines[0].a1/a + splines[0].a0
#b = Math.log(splines[0].a0 - c) + a*splines[0].x0
#p "%20.12E %20.12E %20.12E" % [a, b, c]
#exit

int_dis = x0s[0]*0.5
npoint = 100
int = (x0s[0]-int_dis)/npoint
(0..npoint-1).each do |j|
    r = int_dis+int*j

    val = 0
    if (deriv == 0)
        val = Math.exp(-a*r+b)+c
    elsif (deriv == 1)
        val = -a*Math.exp(-a*r+b)
    elsif (deriv == 2)
        val = a*a*Math.exp(-a*r+b)
    end
    if (val.abs > 5)
        next
    end
    
    puts "%.6f %.6f" % [r, val]
end

x0s.each_with_index do |x0, i|
    if (i==x0s.size-1)
        break
    end
    x1 = x0s[i+1]
    int = (x1-x0)/10.0
    (0..9).each do |j|
        r = x0+int*j
        puts "%.6f %.6f" % [r, splines[i].eval(r, deriv)]
    end
    
end
