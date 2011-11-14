#!/usr/bin/env ruby
require 'matrix'

xyz = ARGV[0]

symbols = []
coords = []
File.open(xyz) do |file|
    nat = file.gets.to_i
    file.gets
    nat.times do 
        line = file.gets
        arr = line.split
        symbols << arr[0].strip
        coords << Vector.elements(arr[1..3].map{|x| x.to_f})
    end 
end

edges = []

(0...coords.size).each do |i|
    (0...i).each do |j|
        if ( symbols[i].downcase != "c" or symbols[j].downcase != "c" ) then
            next
        end
        if ( (coords[i]-coords[j]).r < 1.7 ) then
            edges << [i+1, j+1]  
        end 
    end
end

puts "{ #{edges.map {|item| item.join("<->")}.join(",")} }"

