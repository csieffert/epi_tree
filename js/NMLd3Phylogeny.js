var minCollectionDate;
var maxCollectionDate;

//Function to load the newick file from disk and begin building the phylogenetic tree
function loadNewick(evt){
	var data;
    var file = evt.target.files[0];
 	var reader = new FileReader();
 	reader.onload = function(){
 		buildTree(reader.result);
 	}
 	var text = reader.readAsText(file);
}

//handles the input of meta data to compliment the phylogenetic tree.
function loadMeta(evt) {
	var data;
    var file = evt.target.files[0];
    Papa.parse(file, {
      header: true,
      complete: function(results) {  
        data = results.data;
        $("#provinceLegend").show();
        //add the meta data to the side of the tree
        add_location_and_timeline_column(data);
      }
    });
  } 

//Build the phylogenetic tree using the supplied newick file
function buildTree(treeText){
	var newick = Newick.parse(treeText);
	var newickNodes = []

	function buildNewickNodes(node, callback) {
		newickNodes.push(node)
		if (node.branchset) {
			for (var i=0; i < node.branchset.length; i++) {
				buildNewickNodes(node.branchset[i])
			}
		}
	}
	buildNewickNodes(newick)

	d3.phylogram.build('#phylogram', newick, {
		width: 3000,
		height: 4000,
	});
	alterInnerNodeText();
}

//function that will alter the font color describing total SNP's for inner node branches
function alterInnerNodeText(){
	var innerNodes = document.getElementsByClassName("inner node");
	var x =0;
	while(innerNodes[x]){
		//Center the total SNP text element for the inner nodes on the tree
		innerNodes[x].childNodes[0].setAttribute("dx", innerNodes[x].childNodes[0].getAttribute("dx")-10);
		if(innerNodes[x].childNodes[0].textContent === "0"){
			innerNodes[x].childNodes[0].remove();
		}
		var node = innerNodes[x].childNodes[0];
		$(node).attr('fill', '#000');
		$(node).attr('font-size', '14px');
		x++;
	}
}

//Add the timeline column to the current phylogenetic tree based on metadata provided
function add_location_and_timeline_column(data){
	//grab all of the leaf nodes in the current tree and append the appropriate column describing the nodes sampling province and date
	var leafNodes = document.getElementsByClassName("leaf node");
	var x =0;
	minCollectionDate = findMinDate(data);
	maxCollectionDate = findMaxDate(data);
	
	var minParsed = minCollectionDate.toDateString().split(' ');
	var maxParsed = maxCollectionDate.toDateString().split(' ');
	$("#collectionDate-firstSample").append(minParsed[1]+" "+minParsed[3]);
	$("#collectionDate-lastSample").append(maxParsed[1]+" "+maxParsed[3]);

	var totalDays = daydiff(minCollectionDate, maxCollectionDate);
	while(leafNodes[x]){
		var node = leafNodes[x].childNodes[1];
		var nodeStrain = leafNodes[x].childNodes[1].textContent.split(" ")[0];
		var strainFound = 0;
		//check if the strain is present in the meta data file, and if so extract the values for its location and time
		var y =0 ;
		while(data[y]  && !strainFound){
			var trimmedSubmitted = $.trim(data[y].SubmittedNumber);
			if(("'"+data[y].NLEP+"'" == nodeStrain) || (data[y].NLEP == nodeStrain) || (trimmedSubmitted == nodeStrain)){
				//create the bar graph that depicts sample location and date of collection:
				var newdiv = document.createElement( "div" );
				$(newdiv).css("class", "locationTimeBar");
				$(newdiv).css("position", "absolute");
				if(data[y].Province){
					$(newdiv).css("backgroundColor", determineProvince(data[y].Province));
					$(leafNodes[x].childNodes[1]).attr('fill', determineProvince(data[y].Province));
				}
				$(newdiv).css("opacity", '0.5');
				var top = $(leafNodes[x]).offset().top + 4;
				$(newdiv).css('top', top);
				$(newdiv).css('left', '50%');
				if(data[y].IsolatDate){
					$(newdiv).css('width', calculatePercentageTime(data[y].IsolatDate, minCollectionDate, totalDays));
				}
				$(newdiv).css('height', '8px');
				$(newdiv).css('z-index', '4');
				$(newdiv).append(" ");
				$("#collectionDate").append(newdiv);
				
				//change the node circles to reflect the province of origin
				$(leafNodes[x].childNodes[0]).attr('fill', determinePfgeGroup(data[y]));
				$(leafNodes[x].childNodes[0]).attr('stroke', determinePfgeGroup(data[y]));
				
				strainFound = 1;
			}
			y++;
		}

		//if there is no metadata for the strain, indicate that on the element
		if(!strainFound){
			//create the bar graph that depicts sample location and date of collection:
			var newdiv3 = document.createElement( "div" );
			$(newdiv3).css("class", "locationTimeBar");
			$(newdiv3).css("position", "absolute");
			var top = $(leafNodes[x]).offset().top;
			$(newdiv3).css('top', top);
			$(newdiv3).css('left', '50%');
			$(newdiv3).css('width', '20px');
			$(newdiv3).css('height', '15px');
			$(newdiv3).css('z-index', '4');
			$(newdiv3).css('font-size', '8px');
			$(newdiv3).append("N/A");
			$("#collectionDate").append(newdiv3);
			$(leafNodes[x].childNodes[0]).attr('fill', "grey");
			$(leafNodes[x].childNodes[0]).attr('stroke', "grey");
		}
		
		//add lines that will facilitate lining meta data up with phylo tree
		/*if((x%2) == 0){
			//add div to indicate the origin province
			var newdiv2 = document.createElement( "div" );
			$(newdiv2).css("position", "absolute");
			$(newdiv2).css("backgroundColor", 'black');
			$(newdiv2).css('opacity', '0.1');
			var top2 = $(leafNodes[x]).offset().top;
			top2 = top2 + 3.5;
			$(newdiv2).css('top', top2);
			$(newdiv2).css('left', '0px');
			$(newdiv2).css('width', '100%');
			$(newdiv2).css('height', '4px');
			$(newdiv2).css('z-index', '2');
			$(newdiv2).addClass('contrastLine');
			$(newdiv2).append(" ");
			$("#collectionDate").append(newdiv2);
		}*/
		
		x++;
	}
}


function determineProvince(province){
	switch($.trim(province.toUpperCase())){
		case 'BC':
			return 'green';
		case 'AB':
			return 'orange';
		case 'SK':
			return 'black';
		case 'MB':
			return '#336666';
		case 'ON':
			return 'red';
		case 'QC':
			return 'blue';
		case 'NS':
			return 'purple';
		case 'NB':
			return 'darkolivegreen';
		case 'PE':
			return 'brown';
		case 'NL':
			return 'pink';
		case 'YT':
			return 'goldenrod';
		case 'NWT':
			return 'aqua';
		case 'NU':
			return '#FF3399';
		case 'ONTARIO':
			return 'red';
		default:
			break;
	}
}

function determinePfgeGroup(isolate){
	var outbreak = $.trim(isolate.Outbreak);
	var xbal = $.trim(isolate.PFGEXbaI);
	var blnl = $.trim(isolate.PFGEBlnI);
	if(/*outbreak.contains('1404STH2') &&*/ xbal.contains('STHXAI.0002') && blnl.contains('STHBNI.0002')){
		return 'yellowGreen';
	}
	else if(/*outbreak.contains('1408STH2') &&*/ xbal.contains('STHXAI.0002') && blnl.contains('STHBNI.0015')){
		return 'aqua';
	}
	else if(/*outbreak.contains('1408STH2B18') &&*/ xbal.contains('STHXAI.0002') && blnl.contains('STHBNI.0018')){
		return 'orange';
	}
	else{
		return 'grey';
	}
}

function calculatePercentageTime(dateString, minDate, totalDays){
	var date = standardizeDate(dateString);
	var days = daydiff(minDate, date);
	var percentage = days / totalDays;
	var pixels = percentage * 300;
	return pixels + "px";
}

//function to standardize the date format string from the input data
function standardizeDate(dateString){
	if(dateString){
	if(dateString.split('-').length>1){
		date = new Date(dateString.split('-')[0], dateString.split('-')[1], dateString.split('-')[2]);
		return date;
	}
	else if(dateString.split('/').length>1){
		date = new Date(dateString.split('/')[2], dateString.split('/')[1], dateString.split('/')[0]);
		return date;
	}
	else{
		return null;
	}
	}
}

//finds min date in the sample
function findMinDate(data){
	var x = 0;		
	var date = standardizeDate(data[1].IsolatDate);
	while(data[x]){
		var newDate = standardizeDate(data[x].IsolatDate);
		if(newDate && (newDate < date)){
			date = newDate;
		}
		x++;
	}
	return date;
}
//finds max date in the sample
function findMaxDate(data){
	var x = 0;
	var date = new standardizeDate(data[0].IsolatDate);
	while(data[x]){
		var newDate = new standardizeDate(data[x].IsolatDate);
		if(newDate && (newDate > date)){
			date = newDate;
		}
		x++;
	}
	
	return date;
}

//calculates the number of days between two dates
function daydiff(first, second) {
    return (second-first)/(1000*60*60*24);
}

//remove the contrast lines from the view
function removeContrast(){
	$(".contrastLine").remove();
}

//explanation for the source of the data:
function sourceExplanation(){
	var text = "The provincial source of each isolate is described by the fill color for each leaf node on the tree.  Isolates that lack meta data are colored greenYellow.";
	alert(text);
}

//explanation for the collection date bar:
function collectionDateExplanation(){
	var text = "Each bar represents the time of sampling for the associated isolate on the tree.  The range of dates is as follows:\n" +
			"\t Minimum date: Date of earliest sample from group (smallest bar) => "+minCollectionDate+"\n" +
			"\t Max date: Date of latest sampling from group (largest bar) => "+maxCollectionDate;
	alert(text);
}