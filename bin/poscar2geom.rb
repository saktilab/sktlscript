#!/usr/bin/env ruby


require 'matrix'

lines = []
while line = gets
    lines << line
end  


title = lines[0]

factor = lines[1].to_f
basis = lines[2..4].map{|x| x.split.map{|y| y.to_f}}

typenames = lines[5].split
nats = lines[6].split.map{|x|x.to_i}




mode_str = lines[7]
if mode_str[0..0].upcase == "D"
    mode = "Relative"
else
    mode = "Cart"
end

coords = []
lines = lines[8..-1]
nats.each_with_index do |n, i|
    n.times do 
        coords << [i+1, lines.shift]
    end
end



final = basis.map{|line| (Vector.elements(line)*factor).to_a}



puts "
    TypeNames = { #{typenames.join(" ") } }
    TypesAndCoordinates [#{mode}] = {
"
    coords.each do |line|
        puts "        #{line.join("  ")}"
    end
puts "    }
    Periodic = Yes
    LatticeVectors [Angstrom] = {
"
    final.each do |line|
        print "        "
        puts line.map{|x| "%20.14f" % x}.join("  ")
    end
puts "    }
"


