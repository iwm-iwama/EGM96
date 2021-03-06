#!/usr/bin/env ruby
#coding:binary

$VERSION = "iwm20210606"
# <<History>>
#   iwm20210606
#   iwm20200204

#-------------------------------------------------------------------------------
# Convert
#   ("%02d%02d%010.7f", deg, min, sec)
#     Latitude	Longitude
#       383741.358000	2694644.958000
#      -143716.381200	3050116.010400
#       465227.548400	1022655.424400
# To
#   decimal degrees
#     Latitude	Longitude
#       38.628155	269.779155
#      -14.621217	305.021114
#       46.874319	102.448729
#
# <<Environment>>
#   This script executing needs 'Ruby'.
#       http://www.ruby-lang.org/
#
#-------------------------------------------------------------------------------

Signal.trap(:INT) do
	exit
end

def rtnGeoIBLto10B(
	ddmmss # ddmmss.s...
)
	ddmmss = ddmmss.to_f

	sign = 1

	if ddmmss < 0
		sign = -1
		ddmmss = -ddmmss
	end

	sec = ddmmss % 100
	min = ((ddmmss / 100).to_i) % 100
	deg = (ddmmss / 10000).to_i

	return sign * (deg + (min / 60.0) + (sec / 3600.0)).to_f
end

def main()
	$IFN = ARGV[0]

	if ! $IFN
		puts
		puts "Not Input file."
		puts
		puts "ex1 > ruby #{$0} \"Input file\""
		puts "ex2 > ruby #{$0} \"Input file\" > \"Output file\""
		puts
		exit
	end

	if ! File.exist?($IFN)
		puts
		puts "Not exist, '#{$IFN}'"
		puts
		exit
	end

	File.open($IFN, "rb") do
		|fs|
		fs.each_line do
			|ln|
			lat, lng, remarks = ln.strip.split("\t")
			printf("%.6f\t%.6f", rtnGeoIBLto10B(lat), rtnGeoIBLto10B(lng))
			print (remarks ? "\t#{remarks}\n" : "\n")
		end
	end
end

main()
