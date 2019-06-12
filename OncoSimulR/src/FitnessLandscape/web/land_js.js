var rectsel=-1;
var oldcolor;
var newcolor; 
var oldop;

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




function setValues()
 
 	{
	 	document.forms['land_form'].elements['LAND_THRESHOLD'].value = 1;
 		document.forms['land_form'].elements['LAND_OPTCLEAN'].checked =true;
  			
   		document.forms['land_form'].elements['LAND_LOG'].checked =true;
 		document.forms['land_form'].elements['LAND_DRAWTOEND'].value =-1;
		document.forms['land_form'].elements['LAND_DRAWFROM'].value =-1;
		document.forms['land_form'].elements['LAND_REFERENCE'].value =0;

		document.forms['land_form'].elements['LAND_SCALE_W'].value =80;
		document.forms['land_form'].elements['LAND_SCALE_H'].value =100;

		document.forms['land_form'].elements['LAND_ONLY_MUT'].value =-1;
		document.forms['land_form'].elements['LAND_FLAT'].checked =false;
		document.forms['land_form'].elements['LAND_COMPACT'].checked =false;


		
	}
 function myhide(divID) {
    var item = document.getElementById(divID);
    if (item) {
      item.className='cache';
    }
  }
 function myshow(divID) {
    var item = document.getElementById(divID);
    if (item) {
      item.className='decache';
    }
  }
  function cachedecache(divID) {
    var item = document.getElementById(divID);
    if (item) {
      item.className=(item.className=='cache')?'decache':'cache';
    }
  }

  function decache(divID) {
    var item = document.getElementById(divID);
    if (item) {
      item.className=(item.className=='cache')?'decache':'cache';
    }
  }

  function getFile(file) {
  alert (file);
}

function add_stat(thestring)
{
	v=document.getElementById('LAND_STAT');
	v.innerHTML=thestring;
}
function change_cursor1()
{
  	document.body.className = 'wait';
	land_form.submit();
  }
function change_cursor2()
{
  	document.body.className = 'wait';
	land_form2.submit();
  }

