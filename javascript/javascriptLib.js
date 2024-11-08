function password() {
  var testV = 1;
  var pass1 = prompt('Please Enter the BOREAS Data Password:','');
  while (testV < 3) {
    if (pass1 == null) 
      break;
    if (pass1 == "argo_nogo") {
      location = "http://hitman.gsfc.nasa.gov/";
      break;
      } 
    testV+=1;
    var pass1 = prompt('Access Denied - Password Incorrect, Please Try Again.','');
    }
  if (pass1!="argo_nogo" & testV ==3)
     alert('Sorry, you cannot access the BOREAS CD-ROM simulation site.');
}						

