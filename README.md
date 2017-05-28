# Volumes-of-Polytopes
Thesis 2 - GTU
Two notions of polytope volume were combined in this project, namely the
Euclidean volume and the discrete volume. The Euclidean volume is the
classical notion of volume, while the discrete volume is the number of inte-
gral points in the polytope.

Emiris – Fisikopoulos’ fast algorithms were used to approximate the Eu-
clidean volume(note that exact computation of volume is #P hard), and
were combined with Ehrhart Theory to compute the discrete volume. The
polytope specified by the integer values in the file with extension .ine is read,
and rewritten to a file which is a new dilated polytope. The volume approxi-
mation determines the upper and lower bounds of the new system. Binomial
and monomial calculations are made with dilated polytope and b and H val-
ues are determined. Linear programming creates a file with lp extension to
obtain H vector. The Ehrhart polynomial obtained using the GLPK library
and the discrete volume are found using LattE.

Emiris-Fisikopoulos’ fast volume approximation that uses the CGAL library,
have been used for Euclidean volume. All steps were implemented as afore-
mentioned. The system has been tested with different types of polytopes
such as; cube, birkhoff, crosspolytope.

USER MANUEL
• Install the Dokcer into the OS with sudo apt-get install docker.io or
Instructions can be followed from the docker website for installation
• Download the latest version of docker image from the following link to
installation of the necessary tools for project with docker pull yunusem-
regtu/bitirme command;
– https://hub.docker.com/r/yunusemregtu/bitirme/
• Run the docker image with sudo docker run ’imageName’
• Download LattE integrale bundle from the following link and extract
file to build file https://www.math.ucdavis.edu/ latte/software.php.
Respectively, give the ./configure prefix=/usr/local and make com-
mand on command line
• Clone or Download the python code that calculates the volumes from
the following link: https://github.com/yunusemravci/Volumes-of-Polytopes
• Exctract the directory to project file
• Give the python Volumes.py polytopeName.ine dilateNumber command
on the command line
