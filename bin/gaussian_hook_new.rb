#!/usr/bin/env ruby

require 'yaml'
require 'tmpdir'

class String
    def trim sep=/\s/
        sep_source = sep.is_a?(Regexp) ? sep.source : Regexp.escape(sep)
        pattern = Regexp.new("\\A(#{sep_source})*(.*?)(#{sep_source})*\\z")
        self[pattern, 2]
    end
end
def which(cmd)
  raise ArgumentError.new("Argument not a string: #{cmd.inspect}") unless cmd.is_a?(String)
  return nil if cmd.empty?
  case RbConfig::CONFIG['host_os']
  when /cygwin/
    exts = nil
  when /dos|mswin|^win|mingw|msys/
    pathext = ENV['PATHEXT']
    exts = pathext ? pathext.split(';').select{ |e| e[0] == '.' } : ['.com', '.exe', '.bat']
  else
    exts = nil
  end
  if cmd[File::SEPARATOR] or (File::ALT_SEPARATOR and cmd[File::ALT_SEPARATOR])
    if exts
      ext = File.extname(cmd)
      if not ext.empty? and exts.any?{ |e| e.casecmp(ext).zero? } \
      and File.file?(cmd) and File.executable?(cmd)
        return File.absolute_path(cmd)
      end
      exts.each do |ext|
        exe = "#{cmd}#{ext}"
        return File.absolute_path(exe) if File.file?(exe) and File.executable?(exe)
      end
    else
      return File.absolute_path(cmd) if File.file?(cmd) and File.executable?(cmd)
    end
  else
    paths = ENV['PATH']
    paths = paths ? paths.split(File::PATH_SEPARATOR).select{ |e| File.directory?(e) } : []
    if exts
      ext = File.extname(cmd)
      has_valid_ext = (not ext.empty? and exts.any?{ |e| e.casecmp(ext).zero? })
      paths.unshift('.').each do |path|
        if has_valid_ext
          exe = File.join(path, "#{cmd}")
          return File.absolute_path(exe) if File.file?(exe) and File.executable?(exe)
        end
        exts.each do |ext|
          exe = File.join(path, "#{cmd}#{ext}")
          return File.absolute_path(exe) if File.file?(exe) and File.executable?(exe)
        end
      end
    else
      paths.each do |path|
        exe = File.join(path, cmd)
        return File.absolute_path(exe) if File.file?(exe) and File.executable?(exe)
      end
    end
  end
  nil
end
def copy_template_dftbplus_input(source_path, target_path)
    path = File.expand_path(source_path)
    if (!File.directory?(source_path))
        path = File.expand_path(File.dirname(source_path))
    end
    input_file = File.expand_path("dftb_in.hsd", path)

    if (!File.exists?(input_file) or !File.readable?(input_file))
        puts "The template DFTB+ input file does not exist!!"
        puts "Path: #{input_file}"
        3.times {puts}
        exit 1
    end

    FileUtils.cp(input_file, target_path)
    final_input_path = File.expand_path("dftb_in.hsd", target_path)
    succ = File.exists?(final_input_path)
    if (succ)
        puts "Copying dftb_in.hsd from #{path} to #{target_path} ... done"
    else
        puts "Copying dftb_in.hsd from #{path} to #{target_path} ... failed"
        3.times {puts}
        exit 1
    end
    
    puts "Scanning dependency of input file..."
    input = File.read(final_input_path)
    pat = /<<[<+]\s*(?<filename>.*?)\s*$/im
    files = input.scan(pat)
    files.flatten!
    files.map!{|x| x.trim '"'}
   
    puts "#{files.size} dependency found."
 
    files.each do |file|
        full_path = File.expand_path(file, path)
    
        if (!File.exists?(full_path) or !File.readable?(full_path))
            "-"*80
            puts "WARNING!!! The inclusion files of the DFTB+ input file does not exist!!"
            puts "If this inclusion file is the geometry file, then you dont need to worry."
            puts "Path: #{full_path}"
            "-"*80
            next
        end

        FileUtils.cp(full_path, target_path)
        succ = File.exists?(File.expand_path(file, target_path))
        if (succ)
            puts "Copying #{file} from #{path} to #{target_path} ... done"
        else
            puts "Copying #{file} from #{path} to #{target_path} ... failed"
            3.times {puts}
            exit 1
        end
    end
end
def write_gen(filename, elements, coords)
    File.open(filename, "w") do |file|
        file.puts "#{elements.size} C"
        symbols = elements.sort.uniq
        file.puts "#{symbols.join(" ")}"
        coords.each_index do |ind|
            file.puts "#{ind+1} #{symbols.index(elements[ind])+1} #{coords[ind].map{|x|x*0.529177}.join(" ")} "
        end
    end
end


def write_dftbplusinput(tmp_dir, deriv, charge, elements, coords)
    all_symbol = elements.sort.uniq
    input_file_path = File.expand_path("dftb_in.hsd", tmp_dir)
    puts "Writing geometry"
    write_gen(File.expand_path("in.gen", tmp_dir), elements, coords)
    File.open(input_file_path, "a") do |file|
        file.puts "!Geometry = !GenFormat {"
        file.puts '   <<< "in.gen"'
        file.puts "}"
        
        file.puts "+Hamiltonian = +DFTB{"
        file.puts "    *Charge = #{charge}"
        if (deriv == 2)
            file.puts "    *SCCTolerance = 1.0e-12"
        end
        file.puts "}"

        if (deriv == 2)
            file.puts "!Driver = !SecondDerivatives {}"
        else 
            file.puts "!Driver = {}"
            if (deriv == 0)
                 
            elsif (deriv == 1)
                file.puts "*Analysis = { "
                file.puts "    *CalculateForces = Yes"
                file.puts "}"
            end
        end
        file.puts "*ParserOptions = {"
        file.puts "    *IgnoreUnprocessedNodes = Yes"
        file.puts "}"
    end
end





ENV['OMP_NUM_THREADS'] = ENV['NSLOTS']
ENV['TMPDIR']=ENV['GAUSS_SCRDIR']
job_path = nil
if (job_path == nil)
    job_path = ENV['SGE_O_WORKDIR']
end
if (job_path == nil ) 
    job_path = ENV['PBS_O_WORKDIR']
end

puts
puts
puts "="*80
puts

puts "*"*80
puts "Gaussian 09 External Script for DFTB+"
puts "Powered by Chien-Pin Chou"
puts "*"*80

puts
dftb_program = '/users/solccp/bin/dftb+_1.3_static'
if dftb_program == nil
    puts "Cannot locate DFTB+ program from $PATH"
else
    puts "DFTB program is located at #{dftb_program}"
end

script_path = File.expand_path(File.dirname(__FILE__))
puts "This script is located at #{script_path}"
require "#{script_path}/chemdata"
puts "Current path: #{Dir.pwd}"
pwd = Dir.pwd
puts "The path of submission: #{job_path == nil ? "nil" : job_path}"

looking_path = job_path

if (job_path == nil)
    puts "no path of submission => use current path"
    looking_path = Dir.pwd
end
    
tmp_dir = File.expand_path(Dir.mktmpdir("gau_hook-"), ENV['TMPDIR'])
puts "Temp path for DFTB+ calculation is #{tmp_dir}"


puts "Reading inputs from Gaussian..."
layer = ARGV[0]
inputfile = File.open(ARGV[1])
nAtoms, deriv, charge, spin = inputfile.gets.split.map{|x| x.to_i}
coords = []
elements = []
nAtoms.times do
    an, x, y, z, mm_charge = inputfile.gets.split.map{|x| x.to_f}
    an = an.to_i
    symbol = ChemData.query_symbol(an)
    elements << symbol
    coords << [x,y,z]
end
inputfile.close

outputfile_name = ARGV[2]
msgfile_name = ARGV[3]
fchkfile_name = ARGV[4]
matelfile_name = ARGV[5]

msgfile = File.open(msgfile_name, "w")

puts "Copying template inputs..."
copy_template_dftbplus_input(looking_path, tmp_dir)

puts "Preparing the final DFTB+ input..."
write_dftbplusinput(tmp_dir, deriv, charge, elements, coords)

Dir.chdir(tmp_dir)

if deriv == 0
    puts "Request DFTB Energy"
elsif deriv == 1
    puts "Request DFTB Energy and Gradient"
elsif deriv == 2
    puts "Request DFTB Hessian"
end

puts "Executing the DFTB+ program..."
`#{dftb_program} | tee output`
if `grep ERROR output`.size > 0 then
    msgfile.puts "="*80
    msgfile.puts `grep -A 100 ERROR output`
    msgfile.puts "="*80
    exit(2)
end
puts "Loading results from DFTB+ calculation"

grad=`sed -n '/Total Forces/,/^ $/p' detailed.out | sed -e '1d' -e '$d'`.split("\n").map{|x| x.split.map{|y| -(y.to_f)}}
energy = `grep 'Total energy' detailed.out | awk '{print $3}'`
dipole=`grep 'Dipole moment  :' detailed.out | awk 'NR==1{print $4,$5,$6}'`.split.map{|x|x.to_f}

outputfile = File.open(outputfile_name, "w")

outputfile.puts [energy, dipole].flatten.map{|x| "%20.12E" % x}.join("")
#outputfile.puts [energy, [0.0,0.0,0.0]].flatten.map{|x| "%20.12E" % x}.join("")
grad.each do |line|
    outputfile.puts line.map{|x| "%20.12E" % x}.join("")
end

(3*nAtoms+2).times do 
    outputfile.puts [0.0,0.0,0.0].map{|x| "%20.12E" % x}.join("")
end

if deriv==2 then
    hess = []
    File.open("hessian.out") do |file|
        hess=file.readlines.map{|x| x.split.map{|y|y.to_f}}.flatten
    end
    hessian = hess.each_slice(3*nAtoms).to_a
    tri_hess = []
    (1..3*nAtoms).each do |i|
        (1..i).each do |j|
            tri_hess << hessian[i-1][j-1]
        end
    end
    tri_hess.each_slice(3) do |line|
        outputfile.puts line.map{|x| "%20.12E" % x}.join("")
    end
end

outputfile.close

puts `cat #{outputfile_name}`
msgfile.close

Dir.chdir(pwd)

puts
puts "="*80
puts
puts
puts
exit



