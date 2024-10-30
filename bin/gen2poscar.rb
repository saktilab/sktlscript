#!/usr/bin/env ruby

require 'matrix'


genfile = ARGV[0]

type = "C"

symbols = []
index = []
coords = []
nat = 0

lat_vector = []

list_symbol = []
list_symnat = []

File.open(genfile, "r") do |file|
    line = file.gets
    arr = line.split
    if arr.size > 1
        type = arr[1]
    end
    nat = arr[0].to_i
    symbols = file.gets.split
    (1..nat).each do |i|
        line = file.gets
        arr = line.split
        index << arr[1].to_i
        if (list_symbol.size == 0)
            list_symbol << symbols[arr[1].to_i-1]
            list_symnat << 1
        else
            if (symbols[arr[1].to_i-1] == list_symbol[-1])
                list_symnat[-1] +=1
            else
                list_symbol << symbols[arr[1].to_i-1]
                list_symnat << 1
            end
        end
        coords << Vector.elements(arr[2..4].map{|x| x.to_f})
    end
    if ( type.upcase != "C")
        file.gets
        (1..3).each do |i|
            arr = file.gets.split.map{|x| x.to_f}
            lat_vector << Vector.elements(arr)
        end
    end
end
 

formula = []
symnat = []
symbols.each_with_index do |sym, i|
    symnat << index.count(i+1)
    formula << sym << index.count(i+1)
end 


puts formula.join
puts "1.0"
if ( type.upcase != "C")
    puts "#{lat_vector[0].to_a.map{|x| "%16.12f"%x}.join(" ")}"
    puts "#{lat_vector[1].to_a.map{|x| "%16.12f"%x}.join(" ")}"
    puts "#{lat_vector[2].to_a.map{|x| "%16.12f"%x}.join(" ")}"
end
puts list_symbol.join(" ")
puts list_symnat.join(" ")
#puts symbols.join(" ")
#puts symnat.join(" ")
if type.upcase == "F"
    puts "Direct"
else
    puts "Cart"

end
(1..nat).each do |i|
    puts [coords[i-1].to_a.map{|x| "%16.12f"%x}.join(" "),symbols[index[i-1]-1],].to_a.join("  ")
end
