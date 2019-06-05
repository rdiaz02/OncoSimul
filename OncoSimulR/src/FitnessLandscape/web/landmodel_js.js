
var rectsel=-1;
var oldcolor;
var newcolor; 
var oldop;
var olddegrade=-1;

/*-----------------------------------------------------*/
function displayID(id,color,evt,op)
{
var compact=document.forms['land_form'].elements['LAND_COMPACT'].value;

if(evt.altKey)
 {
	document.forms['land_form'].elements['LAND_DRAWFROM'].value = id;
 newcolor='#339999';
 }
else
if(evt.shiftKey)
 {
	document.forms['land_form'].elements['LAND_DRAWTOEND'].value = id;
 newcolor='#996699';
 }
else
 {
	document.forms['land_form'].elements['LAND_REFERENCE'].value = id;
 newcolor='#000000';
 }
	var clickedElement = document.getElementById(id);
	if (rectsel!=-1)
					{
					var svgElement = document.getElementById(rectsel);
					svgElement.style.stroke ="none";
					svgElement.style.fill=oldcolor;
					svgElement.style['stroke-width']=0;
					}
 			rectsel=id;
 			oldcolor=color;
				clickedElement.style.stroke=newcolor;
				clickedElement.style['stroke-width']=2;
		
 }

/*-----------------------------------------------------*/
function displayIDdegrad(id,evt)
{


if(evt.altKey)
 {
	document.forms['land_form'].elements['LAND_DRAWFROM'].value = id;
 newcolor='#339999';
 }
else
if(evt.shiftKey)
 {
	document.forms['land_form'].elements['LAND_DRAWTOEND'].value = id;
 newcolor='#996699';
 }
else
 {
	document.forms['land_form'].elements['LAND_REFERENCE'].value = id;
 newcolor='#000000';
 }
	var clickedElement = document.getElementById(id);
	if (rectsel!=-1)
					{
					var svgElement = document.getElementById(rectsel);
					svgElement.style.stroke ="none";
					svgElement.style['stroke-width']=0;
					}
 	rectsel=id;

	clickedElement.style.stroke=newcolor;
	clickedElement.style['stroke-width']=2;
		
 }

function checklog(v)
{
	if (v<=0)
	{
	alert("Error: Log cannot be chosen when Fitness values <=0 ");
//	  document.getElementById("LAND_LOG").checked = false;
	document.forms['land_form'].elements['LAND_LOG'].checked = false;
	}
}
/*-----------------------------------------------------*/
function checkform()
{
var ok=1;

      
if (document.forms['land_form'].elements['LAND_NBLOCI'].value<2)
	{
	alert ("You must specify a number of loci >=2");
	ok=0;
	}
else
	{

	a = document.getElementById('LAND_NBLOCI').value;
	s=1;
	for (var i = 0; i < a; i ++)
   		{
   		
      	name="LAND_ALL_"+i;
      	if (document.getElementById(name).value<2)
      		{ok=0;	alert ("locus "+(i+1)+": you must specify a value >= 2");break;}
       	s=s*document.getElementById(name).value;	
      	}
     if (s>200)
      	{alert ("Too many genotypes.\nReduce number of loci or/and allels or run in command line mode (see help section)");
    	ok=0;
    	}
     m=document.getElementById('LAND_MODEL').value; 	
     if (m=="LAND_Full_Model")
     {
     }
    }  	
    
if (ok==1)   
	{ 
	LAND_WINDOW_WIDTHSIZE=window.screen.availWidth;
	LAND_WINDOW_HEIGHTHSIZE=window.screen.availHeight
	document.body.className = 'wait';
	land_form.submit();
	}
  	
}

/*-----------------------------------------------------*/
function resetValues(n)
 
 	{
/* 	sortie="";
 	for (var i=0;i< document.forms['land_form'].elements.length;i++)
 	{
		sortie=sortie +"\nel "+i;
	
		sortie=sortie+" = ";
		sortie=sortie+document.forms['land_form'].elements[i].value;
 	}*/
 
	 	document.forms['land_form'].elements['LAND_THRESHOLD'].value = 1.0;


 		document.forms['land_form'].elements['LAND_OPTCLEAN'].checked =true;
   			
 		document.forms['land_form'].elements['LAND_REFERENCE'].value =0;
 		
 		document.forms['land_form'].elements['LAND_DRAWFROM'].value =-1;
  		document.forms['land_form'].elements['LAND_DRAWTOEND'].value =-1;
		
  		document.forms['land_form'].elements['LAND_LOG'].checked =false;
		
		document.forms['land_form'].elements['LAND_CHAINS'].checked =false;
 		document.forms['land_form'].elements['LAND_SCALE_W'].value =80.0;
		document.forms['land_form'].elements['LAND_SCALE_H'].value =100.0;
		document.forms['land_form'].elements['LAND_ONLY_MUT'].value =-1;
		document.forms['land_form'].elements['LAND_FLAT'].checked =false;
		document.forms['land_form'].elements['LAND_COMPACT'].checked =false;
	
	for (var i=0; i<n; i++)
			{
		
			m='LAND_MASK_'+i;
		
			document.forms['land_form'].elements[m].value=-1;
		
			}
/*	sortie=sortie+"\napres changement\n";
 	for (var i=0;i< document.forms['land_form'].elements.length;i++)
 		{
		sortie=sortie +"\nel "+i;
	
		sortie=sortie+" = ";
		sortie=sortie+document.forms['land_form'].elements[i].value;

		}
 	alert(sortie);*/
	}
	
	
/*-----------------------------------------------------*/
function reset_previous()
{

document.forms['land_form'].elements['LAND_PREVIOUS_FILE'].value="";
}	


/*-----------------------------------------------------*/
  function cachedecache(divID) {
    var item = document.getElementById(divID);
    if (item) {
      item.className=(item.className=='cache')?'decache':'cache';
    }
  }
  
 /*-----------------------------------------------------*/ 
 function myhide(divID) {
    var item = document.getElementById(divID);
    if (item) {
      item.className='cache';
    }
  }

/*-----------------------------------------------------*/
 function myshow(divID) {
    var item = document.getElementById(divID);
    if (item) {
      item.className='decache';
    }
  }

/*-----------------------------------------------------*/
function CreateTable()
{
	a = document.getElementById('LAND_NBLOCI').value;
	m=document.getElementById('LAND_MODEL').value;
    var tablecontents = "";
    tablecontents = "<table >Enter your numbers of alleles";
    tablecontents += "<tr>";
    name="LAND_ALL_0";
    tablecontents += "<td><input type=\"text\" id=\""+name+"\" name=\""+name+"\" size=\"1\" onchange=\"UpdateTable()\" value=\"2\"> </td>";

    for (var i = 1; i < a; i ++)
   		{
      	name="LAND_ALL_"+i;
      	tablecontents += "<td><input type=\"text\" id=\""+name+"\" name=\""+name+"\" size=\"1\" > </td>";
     
  		 }
	tablecontents += "</tr>";
   	tablecontents += "</table>";
	tablecontents += "<p>Same number of alleles<input type= \"checkbox\" id=\"LAND_CHECK_SAMENBR\" name=\"LAND_CHECK_SAMENBR\" checked onclick=\"UpdateTable()\"></p>";
 	tablecontents += "<p>	<input type=\"Button\" value= \"Draw\" onClick=checkform(); />	";
   	document.getElementById("tablespace").innerHTML = tablecontents;
   	
  
   	
   	
   	
	UpdateTable();

}

/*-----------------------------------------------------*/
function UpdateTable()
{

	
	v=document.getElementById('LAND_ALL_0').value;
	b=document.getElementById('LAND_CHECK_SAMENBR');
	n=document.getElementById('LAND_NBLOCI').value;
	
	h=n+","+v;
	if(b.checked==true)
		{
		//alert("true "+n);	
   		for (var i = 1; i < n; i ++)
   			{
   			
   			var d="LAND_ALL_"+i;
   			var thecell=document.getElementById(d);
   			h=h+","+v;
   			thecell.value=v;
   			thecell.style.background = "#DDD";
   			thecell.disabled = true;
   			}
   		}
	else
		{

 		for (var i = 1; i < n; i ++)
   			{
   	
   			var d="LAND_ALL_"+i;
   			var thecell=document.getElementById(d);
   			thecell.disabled = false;
   			thecell.style.background = "#FFF";
   			thecell.value="";
			
   			}
		}
//alert("done "+n );		
//document.getElementById('LAND_GENO_OLD').value=h;		
}

/*-----------------------------------------------------*/
function add_stat(thestring)
{
	v=document.getElementById('LAND_STAT');
	v.innerHTML=thestring;
}

/*-----------------------------------------------------*/
function addCaseValue()
{

h=document.getElementById('LAND_MULTSAME');

if (h.checked==false)
	{
	// if already selected then erase the value 

	document.getElementById('LAND_VVV').innerHTML="";
	}
else
	{
	
	showText="";
	showText="Value= <input type= \"text\" id=\"LAND_VALUE_MULTSAME\" name=\"LAND_VALUE_MULTSAME\" value=\"0\" size=\"3\" > (0 for random)";
	document.getElementById(id='LAND_VVV').innerHTML=showText;
	}
}

function setvalue(field2,field1,value)
{
document.getElementById(field2).value = document.getElementById(field1).value;
}
/*-----------------------------------------------------*/
function affiche_suite()
{
	m=document.getElementById('LAND_MODEL').value;

	document.getElementById("LAND_NBLOCI").innerHTML="";
	document.getElementById(id="mult_same").innerHTML="";
	showText="";
	if (m=="LAND_Multiplicative")
		{
			showText="<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH></TH><TH><TH></TH></TR><TR>";
			showText=showText+"<TR><td>s:</Td><td><input type=\"text\" value=\"0.1\" size=3 id=\"LAND_mu_a\" name=\"LAND_mu_a\"  ></Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_a\" name=\"LAND_sigma_a\"></Td>";
			showText=showText+"<TD>DimRet:</TD><TD><Input type=\"text\" value=\"0\" size=3 id=\"LAND_DIM_RET\" name=\"LAND_DIM_RET\" readonly STYLE=\"color: lightgrey;\"></TD>";
//			showText=showText+"<TD></TD><TD><Input  type= \"checkbox\" id=\"LAND_MODELE_LOG\" name=\"LAND_MODELE_LOG\" \"></TD>";
			showText=showText+"</TR>";
			showText=showText+"</Table>";
		
		}
		else
		if (m=="LAND_UniRandom")
			{
			showText="<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH></TR>";
			showText=showText+"<TR><td>HOC</Td><td>0</Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_h\"  name=\"LAND_sigma_h\"></Td>";
			
			showText=showText+"</TR>";
			showText=showText+"</Table>";
			}
		else
		if (m=="LAND_Ising")
			{
			showText="<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH><th>circular<TH></TR>";
			showText=showText+"<TR><td>Ising</Td><td><input type=\"text\" value=\"0\" size=3 id=\"LAND_mu_i\" name=\"LAND_mu_i\"></Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_i\" name=\"LAND_sigma_i\"></Td>";
			showText=showText+"<TD><input type= \"checkbox\" id=\"LAND_CIRCULAR_ISING\" name=\"LAND_CIRCULAR_ISING\" ></TD>	";

			showText=showText+"</TR>";
			showText=showText+"</Table>";

			}
		else
		if (m=="LAND_RMF")
			{
			showText="<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH><th>Dim Ret</TH></TR>";
			showText=showText+"<TR><td>HOC</Td><TD></Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_h\" name=\"LAND_sigma_h\" id=\"LAND_sigma_h\" ></Td><TD></TD>";
			showText=showText+"</TR>";
			showText=showText+"<TR><td>Mult</Td><td><input type=\"text\" value=\"0\" size=3 id=\"LAND_mu_a\" name=\"LAND_mu_a\" ></Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_a\" name=\"LAND_sigma_a\"></Td>";
				showText=showText+"<TD><Input type=\"text\" value=\"0\" size=3 id=\"LAND_DIM_RET\" name=\"LAND_DIM_RET\"></TD>";
			showText=showText+"</TR>";
			showText=showText+"</Table>";

			}
		else

		if (m=="LAND_EggBox")
			{
			showText="<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH><th>/TH></TR>";
			showText=showText+"<TR><td>EggBox</Td><td><input type=\"text\" value=\"0\" size=3 id=\"LAND_mu_e\" name=\"LAND_mu_e\"></Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_e\" name=\"LAND_sigma_e\"></Td>";
		showText=showText+"<TD></TD>";	
			
			showText=showText+"</TR>";
			showText=showText+"</table>";
			}
		else
		if (m=="LAND_Optimum")
			{
			showText="<table>";
			showText=showText+"<TR><td></Td><TD colspan=2><b>Production</b></TD><TD colspan=2><b>Fitness Funct.</b></TD></TR>";
			showText=showText+"<TR><td></Td><td>mean</Td><td>stdev</Td><td>mean</Td><td>stdev </Td></TR>";
			showText=showText+"<TR><td></Td><td><input type=\"text\" value=\"0\" size=3 id=\"LAND_mu_p\" name=\"LAND_mu_p\"></TD><TD><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_p\" name=\"LAND_sigma_p\"></Td>";
			showText=showText+"<td><input type=\"text\" value=\"0\" size=3 id=\"LAND_mu_o\" name=\"LAND_mu_o\"></Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_o\" name=\"LAND_sigma_o\"></Td>";
			showText=showText+"</TR>";
			showText=showText+"</table>";

			}
		else
		if (m=="LAND_Kauffman_Nk")
			{
			showText=showText+"<table><TR><th></TH><th>loci</TH><th>random</TH></TR>";
			showText=showText+"<TR><td>Kauf</Td><td><input type=\"text\" value=\"\" size=3 id=\"LAND_NK_PARAM\" name=\"LAND_NK_PARAM\"></Td><td><input type= \"checkbox\"  name=\"LAND_RAND_KAUF\" id=\"LAND_RAND_KAUF\" checked></Td>";
			showText=showText+"</TR>";
			showText=showText+"</table>";
			}
		else
		if (m=="LAND_Full_Models")
			{


			showText=showText+"<table><TR><th>Name</TH><th>mean</TH><th>stdev</TH><th>various</TH></TR>";
			
			showText=showText+"<TR><td>HOC</Td><td>0</Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_h\"  name=\"LAND_sigma_h\" ></Td><td></td>";
			showText=showText+"</TR>";
			showText=showText+"<tr bgcolor=\"#DFDFDF\"><td></td> <td></td> <td></td> <td></td><td></td></tr>";
			showText=showText+"<tr><td></td><td></td><td></td><td></td><td></td><td></td></tr>";
			showText=showText+"<TR><td>s:</Td><td><input type=\"text\" value=\"0.1\" size=3 id=\"LAND_mu_a\" name=\"LAND_mu_a\"></Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_a\" name=\"LAND_sigma_a\"></Td>";
			showText=showText+"<TD>DimRet:</TD><TD><Input type=\"text\" value=\"0\" size=3 id=\"LAND_DIM_RET\" name=\"LAND_DIM_RET\" readonly STYLE=\"color: lightgrey;\"></TD>";			
			showText=showText+"</TR>";
	showText=showText+"<tr bgcolor=\"#DFDFDF\"><td></td> <td></td> <td></td> <td></td><td></td></tr>";
			showText=showText+"<tr><td></td><td></td><td></td><td>circular</td></tr><td></td>";
			showText=showText+"<TR><td>Ising</Td><td><input type=\"text\" value=\"0\" size=3 id=\"LAND_mu_i\" name=\"LAND_mu_i\"></Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_i\" name=\"LAND_sigma_i\"></Td>";
			showText=showText+"<TD><input type= \"checkbox\" id=\"LAND_CIRCULAR_ISING\" name=\"LAND_CIRCULAR_ISING\" ></TD>	";
	showText=showText+"<tr bgcolor=\"#DFDFDF\"><td></td> <td></td> <td></td> <td></td><td></td></tr>";
			showText=showText+"<tr><td></td><td></td><td></td><td></td><td></td></tr>";
			showText=showText+"<TR><td>EggBox</Td><td><input type=\"text\" value=\"0\" size=3 id=\"LAND_mu_e\" name=\"LAND_mu_e\"></Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_e\" name=\"LAND_sigma_e\"></Td>";
			showText=showText+"<TD></TD>";	
			showText=showText+"</TR>";
		showText=showText+"<tr bgcolor=\"#DFDFDF\"><td></td> <td></td> <td></td> <td></td><td></td></tr>";	
		
		
			showText=showText+"<TR><td></Td><TD colspan=2><b>Production</b></TD><TD colspan=2><b>Fitness Funct.</b></TD></TR>";
			showText=showText+"<TR><td>Optimum</Td><td>mean</Td><td>stdev</Td><td>mean</Td><td>stdev </Td></TR>";
			showText=showText+"<TR><td></Td><td><input type=\"text\" value=\"0\" size=3 id=\"LAND_mu_p\" name=\"LAND_mu_p\"></TD><TD><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_p\" name=\"LAND_sigma_p\"></Td>";
			showText=showText+"<td><input type=\"text\" value=\"0\" size=3 id=\"LAND_mu_o\" name=\"LAND_mu_o\"></Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_sigma_o\" name=\"LAND_sigma_o\"></Td>";

		
		
		
		
			
			

						showText=showText+"<tr bgcolor=\"#777777\"><td></td> <td></td> <td></td> <td></td><td></td></tr>";	
				showText=showText+"<tr><TD></TD><td><b>fitness</b></Td><td></Td><td></Td><td></Td></TR>";			
		showText=showText+"<tr><TD>Minimum:</TD><td> <Input type=\"text\" value=\"0\" size=3 id=\"LAND_FIX_AMOUNT\" name=\"LAND_FIX_AMOUNT\"></td><td></Td></TR>";
						showText=showText+"<tr bgcolor=\"#777777\"><td></td> <td></td> <td></td> <td></td><td></td></tr>";	
			showText=showText+"<tr><td></td><td><b>loci</B></td><td><b>random</B></td><td></td><td></td></tr>";
			

	
			showText=showText+"<TR><td>Kauf</Td><td><input type=\"text\" value=\"-1\" size=3 id=\"LAND_NK_PARAM\" name=\"LAND_NK_PARAM\"></Td><td><input type= \"checkbox\" id=\"LAND_RAND_KAUF\"  name=\"LAND_RAND_KAUF\" checked></Td>";
			showText=showText+"<td></td></TR>";
		showText=showText+"<tr bgcolor=\"#777777\"><td></td> <td></td> <td></td> <td></td><td></td></tr>";		
			
			showText=showText+"</table>";
			

//			myshow("mixedland");					
			
			}
			document.getElementById("mixedland").innerHTML=showText;	
}
function change_cursor()
{
  	document.body.className = 'wait';
	land_form.submit();
  }
function change_cursor2()
{
  	document.body.className = 'wait';
	land_form2.submit();
  }
