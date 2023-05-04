package BLAST::DrawSequenceObjects;

use strict;

use warnings;
no warnings;

use Data::Dumper;

use GD::Simple;
use Class::Struct;



struct ("R2::DrawSequenceObjects" =>
	{
	 'imageWidth'       => '$',  #Width of the output image
	 'imageHeight'      => '$',  #height of the output image
	 'seqLength'        => '$',  #Length of the sequence
	 'pxSeqLine'        => '$',  #Length in pixels of the line representing the input sequence
	 'xOffsetSeq'       => '$',  #Left offset in pixels to start drawing the sequence line
	 'yOffsetSeq'       => '$',  #Y-offset to start drawing the sequence line
	 'heightRectangle'  => '$',  #The height of rectangles
	 'yOffsetLine'      => '$',  #Distance between the sequence line and the lines representing objects (e.g. TMS)
	 'outFormat'        => '$',  #Format of the image that will be generated
	 'outImgFile'       => '$',  #File name where to put the output image
	 'lines'            => '@',  #Array with the sequence positions of lines to draw below the sequence
	 'rectangles'       => '@',  #Array with the sequence positions and labels of the rectangles that will be drawn
	}
       );


#==========================================================================
#Table of integer numbers to colors (can be expanded)
my %colors = ( 1 => 'red',  2 => 'cyan',   3 => 'lime',  4 => 'yellow',  5 => 'pink',
	       6 => 'magenta',  7 => 'blue',   8 => 'green',   9 => 'gray',   10 => 'orange',
	      11 => 'navy', 12 => 'brown', 13 => 'teal',    14 => 'purple', 15 => 'olive');



#==========================================================================
#Override constructors of class variables to validate and initialize values
#==========================================================================

#--------------------------------------------------------------------------
#Width of the output image

sub imageWidth {

  my ($self, $value) = @_;


  my $default_value = 1030;


  if ($value) {

    #validate
    unless ( $value < 250 ) {
      die "Error: Image width can't be shorter than 250px";
    }

    $self->{'R2::DrawSequenceObjects::imageWidth'} = $value;
  }

  unless ($self->{'R2::DrawSequenceObjects::imageWidth'}) {
    $self->{'R2::DrawSequenceObjects::imageWidth'} = $default_value;
  }

  return $self->{'R2::DrawSequenceObjects::imageWidth'};
}



#--------------------------------------------------------------------------
#height of the output image

sub imageHeight {

  my ($self, $value) = @_;


  my $default_value = 40;


  if ($value) {

    #validate
    unless ( $value < 40 ) {
      die "Error: Image height can't be shorter than 40px";
    }

    $self->{'R2::DrawSequenceObjects::imageHeight'} = $value;
  }

  unless ($self->{'R2::DrawSequenceObjects::imageHeight'}) {
    $self->{'R2::DrawSequenceObjects::imageHeight'} = $default_value;
  }

  return $self->{'R2::DrawSequenceObjects::imageHeight'};
}


#--------------------------------------------------------------------------
#Length in pixels of the line representing the input sequence

sub pxSeqLine {

  my ($self, $value) = @_;


  my $default_value = 1000;


  if ($value) {

    #validate
    unless ( $value < 200 ) {
      die "Error: Sequence line can't be shorter than 200px";
    }

    $self->{'R2::DrawSequenceObjects::pxSeqLine'} = $value;
  }

  unless ($self->{'R2::DrawSequenceObjects::pxSeqLine'}) {
    $self->{'R2::DrawSequenceObjects::pxSeqLine'} = $default_value;
  }

  return $self->{'R2::DrawSequenceObjects::pxSeqLine'};
}



#--------------------------------------------------------------------------
#Left offset in pixels to start drawing the sequence line


sub xOffsetSeq {

  my ($self, $value) = @_;


  my $default_value = 15;


  if ($value) {

    #validate
    unless ( $value < 5 ) {
      die "Error: X-offset for sequence line can't be shorter than 5px";
    }

    $self->{'R2::DrawSequenceObjects::xOffsetSeq'} = $value;
  }

  unless ($self->{'R2::DrawSequenceObjects::xOffsetSeq'}) {
    $self->{'R2::DrawSequenceObjects::xOffsetSeq'} = $default_value;
  }

  return $self->{'R2::DrawSequenceObjects::xOffsetSeq'};
}



#--------------------------------------------------------------------------
#Y-offset to start drawing the sequence line


sub yOffsetSeq {

  my ($self, $value) = @_;


  my $default_value = 25;


  if ($value) {

    #validate
    unless ( $value < 5 ) {
      die "Error: Y-offset for sequence line can't be shorter than 5px";
    }

    $self->{'R2::DrawSequenceObjects::yOffsetSeq'} = $value;
  }

  unless ($self->{'R2::DrawSequenceObjects::yOffsetSeq'}) {
    $self->{'R2::DrawSequenceObjects::yOffsetSeq'} = $default_value;
  }

  return $self->{'R2::DrawSequenceObjects::yOffsetSeq'};
}



#--------------------------------------------------------------------------
#Y-offset to substract from yOffsetSeq to generate the height of rectangles

sub heightRectangle {

  my ($self, $value) = @_;


  my $default_value = 10;


  if ($value) {

    #validate
    unless ( $value < 5 ) {
      die "Error: Height of rectangle can't be less than 5px";
    }

    $self->{'R2::DrawSequenceObjects::heightRectangle'} = $value;
  }

  unless ($self->{'R2::DrawSequenceObjects::heightRectangle'}) {
    $self->{'R2::DrawSequenceObjects::heightRectangle'} = $default_value;
  }

  return $self->{'R2::DrawSequenceObjects::heightRectangle'};
}



#--------------------------------------------------------------------------
#Y-offset to add to yOffsetSeq to draw lines below the sequence. That is,
#the distance between the sequence line and the lines representing objects
#(e.g. TMS)


sub yOffsetLine {

  my ($self, $value) = @_;


  my $default_value = 5;


  if ($value) {

    #validate
    unless ( $value < 3 ) {
      die "Error: Distance between the line and the sequence can't me less than 3px.";
    }

    $self->{'R2::DrawSequenceObjects::yOffsetLine'} = $value;
  }

  unless ($self->{'R2::DrawSequenceObjects::yOffsetLine'}) {
    $self->{'R2::DrawSequenceObjects::yOffsetLine'} = $default_value;
  }

  return $self->{'R2::DrawSequenceObjects::yOffsetLine'};
}



#--------------------------------------------------------------------------
#Format of the image that will be generated

sub outFormat {

  my ($self, $value) = @_;


  my $default_value = 'png';


  if ($value) {

    #validate
    unless ($value =~ /^(png|jpg|pdf|eps|tif)$/ ) {
      die "Error: Acceptable image format are: png, jpg, pdf, eps and tif";
    }

    $self->{'R2::DrawSequenceObjects::outFormat'} = $value;
  }

  unless ($self->{'R2::DrawSequenceObjects::outFormat'}) {
    $self->{'R2::DrawSequenceObjects::outFormat'} = $default_value;
  }

  return $self->{'R2::DrawSequenceObjects::outFormat'};
}




###########################################################################
#
#       Main Functions from This Point
#
###########################################################################

#==========================================================================
#FUnction to plot sequence objects (e.g. TMS, motifs, domains, etc.)
#on a protein sequence. Objexts are represented by lines below the sequence
#(e.g. to represent TMS) and rectangles above the sequence (e.g. to
#represent domains).

sub plot_lines_vs_rectangles {

  my $self = shift;

  die "THe length of the sequence is mandatory" unless ($self->seqLength);


  #Get pixel coordinates for TMS
  my $lines      = $self->get_pixels_for_lines();
  my $rectangles = $self->get_pixels_for_rectangles();

#  print Data::Dumper->Dump([$lines, $rectangles], [qw(*lines *rectangles )]);
#  exit;


  my $imgObj = GD::Simple->new($self->imageWidth, $self->imageHeight);

  #Default color for lines
  $imgObj->fgcolor('black');


  #The line representing the protein sequence
  $imgObj->moveTo($self->xOffsetSeq,  $self->yOffsetSeq);
  $imgObj->lineTo($self->xOffsetSeq + $self->pxSeqLine, $self->yOffsetSeq);



  #Draw Rectangles representing motifs, repeats, etc.
  foreach my $rect (@{ $rectangles }) {

    my $color = $colors{$rect->[0]};
    $imgObj->bgcolor($color);


    #The coordinates of the rectangle
    my $xRectStartCoord = $self->xOffsetSeq + $rect->[1];
    my $xRectEndCoord   = $self->xOffsetSeq + $rect->[2];
    my $yRectStartCoord = $self->yOffsetSeq - $self->heightRectangle;
    my $yRectEndCoord   = $self->yOffsetSeq;


    #Calculate coordinates for label (centered on top of rectangle)
    my $xFontCoord      = int (($xRectStartCoord + $xRectEndCoord) / 2);
    my $yFontCoord      = $yRectStartCoord - 1;

#    print Data::Dumper->Dump([$rect, $color, $xRectStartCoord, $RectEndCoord, $yRectStartCoord, $yRectEndCoord],
#			     [qw(*motif *xRectStartCoord *xRectEndCoord *yRectStartCoord *yRectEndCoord)]);
#    <STDIN>;


    #Draw rectangle
    $imgObj->rectangle($xRectStartCoord, $yRectStartCoord, $xRectEndCoord, $yRectEndCoord);

    #Draw label for rectangle
    $imgObj->moveTo($xFontCoord, $yFontCoord);
    $imgObj->string($rect->[0]);
  }



  #Draw lines representing the TMS
  foreach my $line (@{ $lines }) {

    $imgObj->penSize(3,3);

    my $xlineStartCoord = $self->xOffsetSeq + $line->[0];
    my $xlineEndCoord   = $self->xOffsetSeq + $line->[1];
    my $ylineCoord      = $self->yOffsetSeq + $self->yOffsetLine;

#    print Data::Dumper->Dump([$line, $xlineStartCoord, $xlineEndCoord, $ylineCoord],
#			     [qw(*tms  *xlineStartCoord *xlineEndCoord *ylineCoord)]);
#    <STDIN>;

    $imgObj->moveTo($xlineStartCoord, $ylineCoord);
    $imgObj->lineTo($xlineEndCoord, $ylineCoord);
  }


  #Generate the png image
  open (my $image, ">", $self->outImgFile) || die $!;
  binmode $image;
  print $image $imgObj->png;
  close $image;
}








#==========================================================================
#Get the pixel coordinates of biological objects that will be drawn with
#rectangles (e.g. motifs, repeats, domains, etc)

sub  get_pixels_for_rectangles {

  my $self = shift;


  return undef unless ($self->seqLength && ($self->lines && @{ $self->rectangles }) &&  $self->pxSeqLine);


  my @pixelCoords = ();
  foreach my $rect (@{ $self->rectangles }) {
    my ($color, $left, $right) = @{ $rect };

    my $pxLeft  = int $left  *  $self->pxSeqLine / $self->seqLength;
    my $pxRight = int $right *  $self->pxSeqLine / $self->seqLength;

    push (@pixelCoords, [$color, $pxLeft, $pxRight]);
  }

  return \@pixelCoords;

}





#==========================================================================
#Get the pixel coordinates of biological objects that will be drawn with
#lines (e.g. TMS)

sub  get_pixels_for_lines {

  my $self = shift;


  return undef unless ($self->seqLength && ($self->lines && @{ $self->lines }) &&  $self->pxSeqLine);


  my @pixelCoords = ();
  foreach my $line (@{ $self->lines }) {

    my ($left, $right) = @{ $line };

    my $pxLeft  = int $left  * $self->pxSeqLine / $self->seqLength;
    my $pxRight = int $right * $self->pxSeqLine / $self->seqLength;

    push (@pixelCoords, [$pxLeft, $pxRight]);
  }

  return \@pixelCoords;

}


1;
