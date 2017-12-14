package TCDB::CheckDependencies;


use warnings;
use strict;
use Data::Dumper;

use Class::Struct;


#==========================================================================
# The purpose of this class is to verify that all third party programs
# required by a script in order to run properly are installed and
# accessible.  Programs will be searched in:
#
# a) the invironment table $path via the 'which' command
# b) a list of user supplied directories.
#
# If a dependency is not properly validated, the dependencies triggering
# the error will be reported and the program will abort.
#
# The class has the following variables:
#
# dependencies_list
#    Reference to an array of dependencies to be verified.
#
# user_dirs
#    Reference to an array of additional dirs to search for
#    dependencies. This is necessary when some programs are
#    installed in directories not registered in the environment
#    variable $PATH.
#
# allow_no_dependencies
#    Flag that if true (1) won't generate an error if the list
#    of dependencies is empty. Useful when programs have no
#    other dependencies. Default is undef;
#
# allow_nonexecutables
#    Flag that if true (1) won't generate an error if a dependency
#    is found but the program is not executable. Default is undef.
#
#
#
# The function that verifies dependencies is: checkDependencies()
#
#
# -------
#
#
# GENERAL GUIDE FOR USAGE:
#
# Load class:
# use TCDB::CheckDependencies;
#
#
# Constructor:
#
# my $myObj = new TCDB::CheckDependencies();
#
#
#
# Get help about this class:
#
# $myObj->help;
#
#
#
# Initiallizing variables through constructor
#
# my $myObj = new TCDB::CheckDependencies(
#    dependencies_list => ['dep1', 'dep2', 'dep3'...],
#    user_dirs => ['dir1', 'dir2', 'dir3', ... ],
#    allow_no_dependencies => 1,
#    allow_nonexecutables => 1);
#
#
#
# Initializing variables through accessors:
#
# my $myObj = new TCDB::CheckDependencies();
# $myObj->dependencies_list(['dep1', 'dep2', 'dep3'...]);
# $myObj->user_dirs([dir1', 'dir2', 'dir3', ... ]);
# $myObj->allow_no_dependencies(1);
# $myObj->allow_nonexecutables(1);
#
#
#
# Once dependencies are defined, verify that all dependencies exist in
# in the environment:
#
# $myObj->checkDependencies;
#
#
#
# -------
# Date created:  11/13/2015
# Programmer:    Arturo Medrano
# Institution:   UCSD, WLU
#==========================================================================



struct ("TCDB::CheckDependencies" =>
	{
	 'dependencies_list'     => '@',  #Input programs for which dependencies will be tested.
	 'user_dirs'            => '@',   #Additional dirs that will be searched for dependencies.
	 'allow_no_dependencies' => '$',  #Turn this flag to 1 if the script won't have dependencies.
	 'allow_nonexecutables'  => '$',  #Turn this flag to 1 if dependencies may not be executable.
	}
       );







#==========================================================================
# This function verifies that the list of dependencies are installed
# in the systems $path variable and/or in a set of user supplied
# directories.
#


sub checkDependencies {

    my ($self, $ar_progList) =  @_;

    my $output = {};

    #Check if a list of programs were given as argument to this function.
    #If the list of programs is not detected it is assumed that it was given
    #to the object when it was initialized.
    if ($ar_progList) {
	$self->dependencies_list($ar_progList);
    }


    #Check if the list of dependencies to search is not available
    unless (@{ $self->dependencies_list } || $self->allow_no_dependencies) {
      die "Error: you must specify a list of dependencies for this script.\n\n",
	"Source -->";
    }


    #Verify that every program is installed here
  PROG:foreach my $prog (@{ $self->dependencies_list }) {


      #Check if program is in variable $path
      my $whichOutput = $self->iWhich($prog);

      #Dependency in $path and is executable
      if (! $whichOutput ) {
	  next PROG;
      }


      #Dependency in $path but is NOT executable
      elsif ($whichOutput == 102) {
	  push (@{ $output->{102} }, $prog);
	  next PROG;
      }

	
      #To this point the program is not in variable $path, check if the user
      #supplied additional dirs to search
      foreach my $upath (@{ $self->user_dirs }) {

	  my $uProg = "$upath/$prog";

	  #Dependency is in user-supplied directory and is executable
	  if (-x $uProg) {
	      next PROG;
	  }

	  #Dependency is in user-supplied directory but is NOT executable.
	  elsif (-e $uProg) {
	      push (@{ $output->{102} }, $uProg);
	      next PROG;
	  }
      }


      #Dependency was not detected at all
      push (@{ $output->{101} }, $prog);
    }


    #If dependency errors where found, report them to user and exit;
    if ( %{ $output } ) {

      my $errorMsg = "";

      foreach my $error (keys %{ $output }) {

	if ($error == 101) {
	  $errorMsg .= "Error: Make sure the following program(s) is(are) installed in the system: \n  " .
	    join ("\n  ", @{ $output->{$error} }) . "\n\n";
	}	
	elsif ($error == 102) {
	  unless ($self->allow_nonexecutables) {
	    $errorMsg .= "Error: the following program(s) is(are) not executable: \n  " .
	      join ("\n  ", @{ $output->{$error} }) . "\n\n";
	  }
	}
	else {
	  die "Unknown error code: $error -->";
	}
      }

      die $errorMsg, "Source --> " if ($errorMsg);
    }
}



#==========================================================================
#This function verifies that a dependency is in the environment variable $path.
#Fuction returns an integer, which is a code indication the success or failure in
#finding the dependency. These are the  return codes:
#
# 100: Empty dependency
# 101: dependency not detected
# 102: dependency was found but NOT  executable


sub iWhich {

    my ($self, $dependency) = @_;


    #Dependency can't be empty
    return 100 unless ($dependency);


    my @dirs = split (/:/, $ENV{'PATH'});


  PATH:foreach my $dir (@dirs) {

      my $path = "$dir/$dependency";

      #Found and executable
      if (-x $path) {
	  return 0;
      }

      #Found but not executable
      elsif (-e $path) {
	  return 102;
      }
  }

    #Dependency not found
    return 101;
}




#==========================================================================
#The help for this class

sub help {

  my $self = shift;

  my $help = << 'HELP';

==========================================================================
 The purpose of this class is to verify that all third party programs
 required by a script in order to run properly are installed and
 accessible.  Programs will be searched in:

 a) the invironment table $path via the 'which' command
 b) a list of user supplied directories.

 If a dependency is not properly validated, the dependencies triggering
 the error will be reported and the program will abort.

 The class has the following variables:

 dependencies_list
    Reference to an array of dependencies to be verified.

 user_dirs
    Reference to an array of additional directories to search for
    dependencies. This is necessary when some programs are
    installed in directories not registered in the environment
    variable $PATH.

 allow_no_dependencies
    Flag that if true (1) won't generate an error if the list
    of dependencies is empty. Useful when programs have no
    other dependencies. Default is undef;

 allow_nonexecutables
    Flag that if true (1) won't generate an error if a dependency
    is found but the program is not executable. Default is undef.



 The function that verifies dependencies is: checkDependencies()


 -------


 GENERAL GUIDE FOR USAGE:

 Load class:
 use TCDB ::CheckDependencies;


 #Object Constructor:
 my $myObj = new TCDB::CheckDependencies();



 #Get help about this class:
 my $myObj = new TCDB::CheckDependencies()->help;


 $myObj->help;



 #Initiallizing variables within the constructor:
 my $myObj = new TCDB::CheckDependencies(
    dependencies_list => ['dep1', 'dep2', 'dep3'...],
    user_dirs => ['dir1', 'dir2', 'dir3', ... ],
    allow_no_dependencies => 1,
    allow_nonexecutables => 1);



 #Initializing variables through accessors:
 my $myObj = new TCDB::CheckDependencies();
 $myObj->dependencies_list(['dep1', 'dep2', 'dep3'...]);
 $myObj->user_dirs(['dir1', 'dir2', 'dir3', ... ]);
 $myObj->allow_no_dependencies(1);
 $myObj->allow_nonexecutables(1);



 #Once dependencies are defined, verify that all dependencies exist in
 #in the environment:
 $myObj->checkDependencies;



 -------
 Date created:  11/13/2015
 Programmer:    Arturo Medrano
 Institution:   UCSD, WLU
==========================================================================

HELP

  print $help;
  exit;
}




1;

