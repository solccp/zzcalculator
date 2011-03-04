#!/usr/bin/env ruby1.9

require 'bigdecimal'
lines = []
File.open("results") do |file|
  while line0 = file.gets
	lines << line0.chomp.strip
  end
end
while lines.size > 0

    line = lines.shift.split
    if line[0][0] == '('
      demo = BigDecimal.new(line[-1])
      name = line[0]
      while lines.size > 0
	 line = lines.shift
         line2 = line.split
         if line2[0].include?("ring") then
             num = BigDecimal.new(line2[-1])
	     res = (num / demo)*BigDecimal.new("100")
             puts "#{"%-20s" % name} #{"%8s" % line2[0]} #{"%15.6f" % res }"
         else
		lines.unshift(line)
		break
         end
      end
    end

end
