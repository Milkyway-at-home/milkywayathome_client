#!/usr/bin/perl

# Takes in 7 arguments: x and y beginning, x and y end, x and y binsize, filename

$xbeginning = $ARGV[0];
$xend = $ARGV[1];
$xbinsize = $ARGV[2];
$ybeginning = $ARGV[3];
$yend = $ARGV[4];
$ybinsize = $ARGV[5];
$filename = $ARGV[6];

open(DAT,$filename);
@values=<DAT>;
close(DAT);

# Record the max index value
$xmaxindex = 0;
$ymaxindex = 0;

foreach $value (@values) {
        chop($value);
	($xvalue,$yvalue) = split(/\s+/,$value);

        # Calculate the index for the histogram data array
        $xindex = int($xvalue/$xbinsize) - int($xbeginning/$xbinsize);
        if($xindex > $xmaxindex) {
                # Record the new max index value
                $xmaxindex = $xindex;
        }
        $yindex = int($yvalue/$ybinsize) - int($ybeginning/$ybinsize);
        if($yindex > $ymaxindex) {
                # Record the new max index value
                $ymaxindex = $yindex;
        }
        # Add 1 to the appropriate array value
        $histodata[$xindex][$yindex]++;
}

# Loop through and print out all of the data
$xindex = 0;
$yindex = 0;

for($i=$xbeginning;$i<=$xend;$i+=$xbinsize) {
	for($j=$ybeginning;$j<=$yend;$j+=$ybinsize) {
        	# The indexes should correspond
		$number = $histodata[$xindex][$yindex];

		if($number eq '') { $number = "0";}

	        $newi = sprintf("%.3f", $i);
	        $newj = sprintf("%.3f", $j);

     		print $newi . " " . $newj . " " . $number . "\n";
        	$yindex++;
	}
	$yindex = 0;
       	$xindex++;
	print "\n";
}
