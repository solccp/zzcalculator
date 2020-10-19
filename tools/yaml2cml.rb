#!/usr/bin/env ruby1.9.1

require 'builder'
require 'yaml'
require 'matrix'

xyz_file = ARGV[0]
yaml_file = ARGV[1]


def write_cml(str, counters, symbol_map, bond_map)
    nrings = 0
    nrings = str["rings"].size if str["rings"]
 
    db = str["double_bonds"]
    if (db) then
        db.each do |a1, a2|
            bond_map[ [a1,a2].sort ] = 2
        end
    end
    rb = str["rings"]
    if (rb) then
        rb.each do |ring|
            ring.each_cons(2) do |a1, a2|
                bond_map[ [a1,a2].sort ] = "A"
            end
            bond_map[ [ring[0], ring[-1]].sort ] = "A"
        end
    end

 
    File.open("ring#{nrings}.str#{counters[nrings]}.cml", "w") do |file|
        xml = Builder::XmlMarkup.new(target: file, indent: 2)
        xml.molecule(xmlns: "http://www.xml-cml.org/schema") do
            xml.atomArray do
                symbol_map.each do |k,v|
                    xml.atom(id: "a#{k}", elementType: "#{v.capitalize}")
                end
            end
            xml.bondArray do
                bond_map.each do |k,v|
                    xml.bond(atomRefs2: "a#{k[0]} a#{k[1]}", order: "#{v}")
                end
            end
        end
    end
    counters[nrings] += 1
end

bond_map = {}
symbol_map = {}

File.open(xyz_file) do |file|
    coord = []
    symbols = []
    nat = file.gets.to_i
    title = file.gets
    nat.times do
        arr = file.gets.split
        symbol, xyz = arr[0], arr[1..3].map{|x| x.to_f}
        coord << xyz
        symbols << symbol
    end

    symbols.each_index do |i|
        symbol_map[i+1] = symbols[i]
    end

    (1..nat).to_a.combination(2) do |arr|
        if ( Vector.elements(coord[arr[0]-1]) - Vector.elements(coord[arr[1]-1])).r < 1.7 then
            bond_map[ [arr[0], arr[1]].sort ] = 1
        end
    end
end

counters = {}
counters.default = 1

File.open(yaml_file) do |file|
    YAML.each_document(file) do |ydoc|
        if (ydoc["number_of_atoms"] != symbol_map.values.select{|x| x.downcase == "c"}.size) then
            puts "different num of atoms" 
            exit(1)
        end
        write_cml(ydoc, counters, symbol_map, bond_map.dup)
    end
end



