# Run through the cvs status reports and find version differences
# after a merge

BEGIN { exists = 1}

/File:/ {
  if (!exists && wRev > 0 && wRev != "New"){
    print "File ", dir, " did not exist on trunk at time of branch.";
  }
  exists=0;
  wRev=0;  
  dRev=0;
  filename = $2;
}

/Working revision:/ {wRev = $3}

/Repository revision:/ {dir = $4}

/draco_121198/ {
  dRev = substr($3,0,length($3)-1);
  exists = 1;
  if (wRev > 0 && wRev != dRev) {
    print "** File ", filename, " has been revised!";
    print "   In directory ", dir;
    print "   Trunk Revision : " wRev;
    print "   Branch Revision:  " dRev;
  }
}

END  {
  if (!exists && wRev > 0 && wRev != "New"){
    print "File ", dir, " did not exist on trunk at time of branch.";
  }
} 
