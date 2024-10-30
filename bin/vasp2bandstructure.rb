#!/usr/bin/env ruby

#data = `sed -n -e '/<set comment="kpoint .*">/,/set/p' vasprun.xml`

#data.each do |x|
#    if x.include?("set comment") 
#    elsif x.include?("set")
#        next
#    end
    
#end

efermi_file = ARGV[0]



require 'rexml/document'
data = `cat vasprun.xml`
fermi = 0.0
if efermi_file!= nil
    fermi = Float(`grep 'E-fermi' #{efermi_file} | tail -n 1 | awk '{print $3}'`)
else
    fermi = Float(`grep 'E-fermi' OUTCAR | tail -n 1 | awk '{print $3}'`)
end
doc = REXML::Document.new(data)
spins = []
alpha_e = 0.0
beta_e = 0.0
doc.elements.each('modeling/calculation/eigenvalues/array/set/set') do |spin|
#should be in spin tag
    spin_index = Integer(spin.attributes["comment"].split[1])
    bands = []
    spin.elements.each('set') do |ele|    
        index = Integer(ele.attributes["comment"].split[1])
        band = []
        ele.elements.each do |r|
            occ = Float(r.text.split[1])
            if (spin_index ==1 )
                alpha_e += occ
            else
                beta_e += occ
            end
            band << Float(r.text.split[0])-fermi
        end
        bands << ["%4d" % index,band.flatten.map{|x| "%10.6f" % x}]
    end
    spins << bands
end
spins.each_with_index do |spin, i|
    File.open("bands_#{i+1}.dat", "w") do |file|
        spin.each do |da|
            file.puts da.join("  ")
        end
    end
end
puts alpha_e, beta_e
