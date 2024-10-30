#!/usr/bin/env ruby

require 'bigdecimal'


class IntergralInfo
    INDEXES     = [ ["d-d-sigma", 1], ["d-d-pi", 2], ["d-d-delta", 3],
                    ["p-d-sigma", 4], ["p-d-pi", 5],
                    ["p-p-sigma", 6], ["p-p-pi", 7],
                    ["s-d-sigma", 8],
                    ["s-p-sigma", 9],
                    ["s-s-sigma",10] ]

    INDEXES_EXT = [ ["f-f-sigma", 1], ["f-f-pi", 2], ["f-f-delta", 3], ["f-f-phi", 4],
                    ["d-f-sigma", 5], ["d-f-pi", 6], ["d-f-delta", 7],
                    ["d-d-sigma", 8], ["d-d-pi", 9], ["d-d-delta",10],
                    ["p-f-sigma",11], ["p-f-pi",12],
                    ["p-d-sigma",13], ["p-d-pi",14],
                    ["p-p-sigma",15], ["p-p-pi",16],
                    ["s-f-sigma",17],
                    ["s-d-sigma",18],
                    ["s-p-sigma",19],
                    ["s-s-sigma",20] ]

    def IntergralInfo.index2name(index, ext_format = false)
        if index < 1
            raise "Index must be > 0."
        elsif index > 10 and !ext_format
            raise "Index of d-integral must be <=10."
        elsif index > 20 and ext_format
            raise "Index of f-integral must be <=20."
        else
            if ext_format
                res = INDEXES_EXT.rassoc(index)[0]
            else
                res = INDEXES.rassoc(index)[0]
            end
        end
        return res
    end

    def IntergralInfo.name2index(name, ext_format = false)
        if ext_format
            res = INDEXES_EXT.assoc(name)
        else
            res = INDEXES.assoc(name)
        end
        return 0 if res == nil
        return res[1]
    end

    def IntergralInfo.filename2index(filename, ext_format = false)
        str = filename.split('-')
        index = name2index([str[5],str[9],str[1]].join("-"), ext_format)
        if ext_format
            index += INDEXES_EXT.size if str[0] == "overl" and index > 0
        else
            index += INDEXES.size if str[0] == "overl" and index > 0
        end
        index
    end
end










def parseLine(line, expand = true)
    res = []
    tmp = line.split(/[,\s]/)
    if (expand == false)
        return tmp
    end
    tmp.each do |x|
        if x==""
            next
        end
        if x.include?("*")
            arr = x.split("*")
            if (arr.size == 2)
                arr[0].to_i.times do
                    res << arr[1]
                end
            end
        else 
            res << x
        end
    end
    res
end


filename = ARGV[0]


integrals = {}
skipping_index = []
File.open(filename, "r") do |file|
    line_1 = parseLine(file.gets.chomp)
    interval = BigDecimal.new(line_1[0])
    zero = BigDecimal.new("0.0")
    line_2 = parseLine(file.gets.chomp)
    if line_2.size >= 10  #onsite
        parseLine(file.gets.chomp)
    end

    dis = interval   
    start_dis = BigDecimal.new("0.4")
    first_line = true
    while line = file.gets
        if line.include?("Spline")
            break
        end
        entries = parseLine(line.chomp, false)
        if entries.size < 4
            dis += interval
            next
        end
        if (dis >= start_dis) 
            entries = parseLine(line.chomp)
            if (first_line)
                first_line = false
                entries.each_with_index do |item, i|
                    if BigDecimal.new(item) == zero
                        skipping_index << i
                    else
                        integrals[i] = Array.new
                        integrals[i] << [dis.to_f, item.to_f]
                    end
                end
            else
                entries.each_with_index do |item, i|
                    if (skipping_index.include?(i))
                        next
                    else
                        integrals[i] << [dis.to_f, item.to_f]
                    end
                end
            end
        end
        dis += interval
    end 
end

coeffs = {}

#integrals.each do |key, value|
#    if (key>9)
#        if (value[0][1]<0.0)
#            coeffs[key-10] = -1.0
#        else
#            coeffs[key-10] = 1.0
#        end
#    end
#end



integrals.each do |key, value|
    index = key+1
    type = "hamil"
    if index > 10
        index -= 10
        type = "overl"
    end
    coeff = 1.0 #coeffs[index-1]
    iname = IntergralInfo.index2name(index)
    prefix = File.basename(filename)
    File.open("#{prefix}.#{type}.#{iname}", "w") do |file|
        value.each do |line|
            file.puts "#{line[0]}  #{line[1]*coeff}"
        end
    end

end
