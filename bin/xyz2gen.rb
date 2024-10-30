#!/usr/bin/env ruby

input_file = ARGV[0]
output_file = `basename #{input_file} .xyz`.chomp.strip + ".gen"

lines = []
symbols = []
natoms = 0
File.open(input_file) do |infile|
    natoms = infile.gets.to_i
    title = infile.gets
    natoms.times do
        line = infile.gets
        arr = line.split
        symbols << arr[0]
        lines << arr
    end 
end

symbols.sort!.uniq!

File.open(output_file, "w") do |outfile|

    outfile.puts [natoms, "C"].join(" ")
    outfile.puts symbols.join(" ")
    lines.each_with_index do |line, index|
        outfile.puts [index+1, symbols.index(line[0])+1, line[1..3]].join("  ")
    end
end
