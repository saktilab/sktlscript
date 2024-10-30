#!/usr/bin/env ruby


require 'matrix'

lines = []
while line = gets
    lines << line
end  


title = lines.shift

factor = lines.shift.to_f
basis_vectors = []
3.times do
   basis_vectors << lines.shift
end

basis = basis_vectors.map{|x| x.split.map{|y| y.to_f*factor}}

typenames = lines.shift.split
nats = lines.shift.split.map{|x|x.to_i}


selective_dynamics = false

if (lines[0].downcase().start_with?("selective"))
    selective_dynamics = true
    lines.shift
end

mode_str = lines.shift
mode_type = "S"
if mode_str[0..0].upcase == "D"
    mode_name = "Relative"
    mode_type = "F"
else
    mode_name = "Angstrom"
end

coords = []
nats.each_with_index do |n, i|
    n.times do 
        line = lines.shift
        if mode_type == "F"
            frac_coord = line.split[0..2].map{|x| x.to_f}
            vec_basis = basis.map{|x| Vector.elements(x)}
            coord = vec_basis[0]*frac_coord[0] + vec_basis[1]*frac_coord[1] + vec_basis[2]*frac_coord[2]
            coords << [typenames[i], coord.to_a.map{|x| "%20.12f" % x}.join("  ")]
        else
            coords << [typenames[i], line.split[0..2].map{|x| "%20.12f" % x}.join(" ")]
        end
    end
end


final = basis.map{|line| (Vector.elements(line)).to_a}

nat = nats.inject(0){|sum,x| sum + x }

puts "#{nat}"
puts ""
coords.each_with_index do |line, j|
    puts "#{line.join("  ")}"
end

#final.each do |line|
#    puts "TV #{line.map{|x| "%20.12f" % x}.join("  ")}"
#end
#puts ""
