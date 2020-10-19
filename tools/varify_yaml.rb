#!/usr/bin/env ruby1.9.1

require 'yaml'

def process(ydoc)

    dbs = []
    rs = []
    str = {}
    str["dbs"] = []
    str["rings"] = []
    if ydoc["double_bonds"] then
        ydoc["double_bonds"].each do |item|
            dbs << item.sort
        end
        str["dbs"] = dbs.sort
    end 
    if ydoc["rings"] then
        ydoc["rings"].each do |item|
            rs << item.sort
        end
        str["rings"] = rs.sort
    end 
    str
end




new_file = ARGV[0]
ref_file = ARGV[1]

ref_strs = {}

File.open(ref_file) do |file|
    YAML.each_document(file) do |ydoc|
        a = process(ydoc)
        ref_strs[a] = 1
    end
end
new_strs = []

File.open(new_file) do |file|
    YAML.each_document(file) do |ydoc|
        a = process(ydoc)
        new_strs << a
    end
end

not_found = []
new_strs.each do |str|
    if ref_strs[str] != 1 then
        not_found << str
    else
        ref_strs.delete(str)
    end
end

p not_found.size
p ref_strs.size

p not_found if not_found.size > 0
p ref_strs if ref_strs.size > 0

