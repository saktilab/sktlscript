#!/usr/bin/env ruby


require 'matrix'


valence_elecs = { :Pt => 8.0, :C => 2.0, :N => 2.5, :O => 3.0, :H => 0.5, :Br => 3.5 }


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
mode_type = "S"
if mode_str[0..0].upcase == "D"
    mode_name = "Relative"
    mode_type = "F"
else
    mode_name = "Angstrom"
end

coords = []
lines = lines[8..-1]
nats.each_with_index do |n, i|
    n.times do 
        line = lines.shift
        if mode_type == "F"
            frac_coord = line.split[0..2].map{|x| x.to_f}
            vec_basis = basis.map{|x| Vector.elements(x)}
            coord = vec_basis[0]*frac_coord[0] + vec_basis[1]*frac_coord[1] + vec_basis[2]*frac_coord[2]
            coords << [typenames[i], coord.to_a.map{|x| "%20.12f" % x}.join("  ")]
        else
            coords << [typenames[i], line.split[0..2].join("    ")]
        end
    end
end


final = basis.map{|line| (Vector.elements(line)*factor).to_a}

nat = nats.inject(0){|sum,x| sum + x }

puts "Atoms.Number #{nat}"
puts "Atoms.SpeciesAndCoordinates.Unit   Ang"
puts "<Atoms.SpeciesAndCoordinates"
coords.each_with_index do |line, j|
    puts "#{j+1}  #{line.join("  ")} #{valence_elecs[line[0].to_sym]} #{valence_elecs[line[0].to_sym]}"
end
puts "Atoms.SpeciesAndCoordinates>"
puts "Atoms.UnitVectors.Unit             Ang"
puts "<Atoms.UnitVectors"
final.each do |line|
    puts "#{line.map{|x| "%20.12f" % x}.join("  ")}"
end
puts "Atoms.UnitVectors>"
puts ""
