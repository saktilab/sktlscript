#!/usr/bin/ruby

file    = ARGV[0]
from_to = ARGV[1] # same phrase to cut out

print_flg = false

File.foreach(file){|line|

	if line =~ /#{from_to}/ and print_flg == false then
		print_flg = true
		next
	end

	if line =~ /#{from_to}/ then
		print_flg = false
		next
	end

	if print_flg then
		print line
	end
}

