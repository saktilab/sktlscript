#!/usr/bin/env ruby

arr = []


filename = ARGV[0]
band_index = ARGV[1].to_i
point_index = ARGV[2].to_i

File.open(filename) do |file|
    while line = file.gets
        a = line.split
        arr << [a[0].to_i, a[1..-1].map{|x| x.to_f}].flatten
    end
end
arr = arr.transpose

shift = arr[band_index+1][point_index]
arr.each_index do |col|
    next if col == 0
    arr[col] = arr[col].map{|x| x-shift}
end

arr = arr.transpose
arr.each do |row|
    puts ["%-5d" % row[0],row[1..-1].map{|x| "%8.4f" % x}].join("  ")
end
