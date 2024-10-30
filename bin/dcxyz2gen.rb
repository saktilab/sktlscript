#!/usr/bin/env ruby

input_file = ARGV[0]
output_file = `basename #{input_file} .xyz`.chomp.strip + ".gen"

lines = []
symbols = []
lattices = []
natoms = 0
File.open(input_file) do |infile|
    natoms = infile.gets.split[0].to_i
    natoms.times do
        line = infile.gets
        arr = line.split
        symbols << arr[0]
        lines << arr
    end
    3.times do 
        line = infile.gets
        arr = line.split
        if (arr[0].downcase == 'tv')
            lattices << arr[1..3]
        end
    end
end

symbols.sort!.uniq!

File.open(output_file, "w") do |outfile|

    if lattices.size == 3
        outfile.puts [natoms, "S"].join(" ")
    else
        outfile.puts [natoms, "C"].join(" ")
    end
    outfile.puts symbols.join(" ")
    lines.each_with_index do |line, index|
        outfile.puts [index+1, symbols.index(line[0])+1, line[1..3]].join("  ")
    end
    if lattices.size == 3
        outfile.puts "0.0 0.0 0.0"
        lattices.each do |line|
            outfile.puts [line[0..2]].join("  ")
        end
    end
end
