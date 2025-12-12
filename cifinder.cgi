#!/usr/local/bin/perl
use strict;

#*******************************************************
# Copyright @2025, Alexei Sharov
# All rights reserved.
# 
# This software is provided "AS IS". Alexei Sharov (AS) makes no warranties, express
# or implied, including no representation or warranty with respect to
# the performance of the software and derivatives or their safety,
# effectiveness, or commercial viability.  AS does not warrant the
# merchantability or fitness of the software and derivatives for any
# particular purpose, or that they may be exploited without infringing
# the copyrights, patent rights or property rights of others. NIA shall
# not be liable for any claim, demand or action for any loss, harm,
# illness or other damage or injury arising from access to or use of the
# software or associated information, including without limitation any
# direct, indirect, incidental, exemplary, special or consequential
# damages.
# 
# This software program may not be sold, leased, transferred, exported
# or otherwise disclaimed to anyone, in whole or in part, without the
# prior written consent of AS.
#
# Programmer: Alexei Sharov (sharov@comcast.net)
#*******************************************************

my %hashIni;
read_configuration("../../cisfinder.ini",\%hashIni);
my $PATH_HOME=$hashIni{"PATH_HOME"};
my $PATH_PROG=$hashIni{"PATH_PROG"};
my $CGI_ADDRESS=$hashIni{"CGI_ADDRESS"};
my $HOME_ADDRESS=$hashIni{"HOME_ADDRESS"};
my $SECURITY=$hashIni{"SECURITY"};
my $UNIX=$hashIni{"UNIX"};
my $PATH_INFO = "$PATH_PROG/info";
my $PATH_DATA = "$PATH_PROG/data";
my $PATH_OUTPUT = "$PATH_HOME/output";
my $LINE = 2;
my $TEXT = 10;
my $POINT = 14;
my $CIRCLE = 1;
my $FILL = 8;
my $radius = 2;
my $black = 15;
my $gray = 14;
my $blue = 8;
my $magenta = 6;
my $red = 21;
my $green = 16;
my $ltgreen = 31;
my $tinyFont = 0;
my $smallFont = 1;
my $largeFont = 2;
my $smallFontInbox = 4;
my $largeFontInbox = 5;
my $vertTinyFont = 6;
my $vertSmallFont = 7;
my $vertLargeFont = 8;
my $printNumbers = 1;
my $thin = 1;
my $thick2 = 2;
my $solid = 0;
my @headersGlob = ("Pattern","Threshold","Nmembers","Freq","Ratio","Info","Score","FDR","Repeat","Palindrome","Method","Species");
my @codes=("O","A","C","M","G","R","S","V","T","W","Y","H","K","D","B","N");
my $CHROM_LENGTH = 50000000;
my $BUFFER = 3000000;

# Read data from the input form
my %hashInput;
my $request_method = $ENV{'REQUEST_METHOD'};
my $size = $ENV{'CONTENT_LENGTH'};
my $form_info;
if($size > 0){
	read (STDIN, $form_info, $size);
} else {
	$form_info = $ENV{'QUERY_STRING'};
}
if(!$form_info){
	$form_info = $ARGV[0];
}
my $loginname;
my $passwd;
my $sessionID;
if(!$form_info){ print "content-type: text/html","\n\n"; print "ERR: no data\n"; exit(0); }
my $date = `date \'+%y%m%d\'`;
chop($date);
if($form_info =~ /\n/){
	print "content-type: text/html","\n\n";
	parse_data(\$form_info,\%hashInput);
	foreach my $key (keys %hashInput){
		if($key=~/</ || $hashInput{$key}=~/</){ error_message("Hyperlinks are blocked!","continue"); }
	}
	$loginname = $hashInput{"loginname"};
	$sessionID = $hashInput{"sessionID"};
	check_sessionID($loginname,$sessionID);
	my $file_type = $hashInput{"file_type"};
	my $rename = $hashInput{"rename"};
	my ($filename,$list);
	my $ref = $hashInput{"upload_file_name"};
	if($ref){
		($filename,$list) = @{$hashInput{"upload_file_name"}};
	}
	if($rename){
		$filename = $rename;
	}
	if($filename =~ /</){ error_message("Hyperlinks are blocked!","continue"); }
	$filename =~ s/.+\///;
	$filename =~ s/.+\\//;
	$filename =~ s/\.(fa|txt)$//;
	$filename =~ s/\./_/g;
	if(!$filename){ error_message("Filename not accepted!","continue"); }
	if($file_type eq "fasta"){ $filename .= ".fa"; }
	elsif($file_type =~ /motif|pattern/){ $filename .= "_motif"; }
	elsif($file_type =~ /search/){ $filename .= "_search"; }
	elsif($file_type =~ /repeat/){ $filename .= "_repeat"; }
	if($file_type =~ /fasta|subset|subregion/ && $filename !~ /\.fa$/){
		$filename .= ".fa";
	}
	my $list1 = $hashInput{"pasted_text"};
	if($list1){
		$list = $list1;
	}
	if($list =~ /</){ error_message("Hyperlinks are blocked!","continue"); }
	my $description = $hashInput{"description"};
	my @lines = split(/\n/,$list);
	if($file_type eq "pattern"){
		if(!pattern_to_PFM(\@lines)){ error_message("Wrong file format!"); }
		$file_type = "motif";
	}elsif($file_type !~ /subregion/){
		if(!check_upload(\@lines,$file_type)){ error_message("Wrong file format!"); }
	}
	$list = join("\n",@lines);
	my $fasta_file = $hashInput{"fasta_file"};
	if($file_type eq "coordinates" && !$fasta_file){
		if($hashInput{"mainPage"} && $hashInput{"action"} =~ /continue|motif_get|motif_search|motif_improve|motif_compare|motif_cluster/){
			check_configuration();
		}
		my ($genome,$junk);
		my $line = $lines[0];
		my @items = split(/\t/,$line);
		while(@items && $items[0] !~ /genome/){ splice(@items,0,1); }
		if($items[0] !~ /genome/){
			error_message("The genome need to be specified in the first line\n(e.g., 'genome=mm9' or 'genome=hg19')");
		}
		($junk,$genome) = split(/=/,$items[0]);
		if($genome ne "mm9" && $genome ne "hg19" && $genome ne "hg18"){
			error_message("Sorry, CisFinder extracts sequences only from mouse genome mm9, or human hg19, hg18.\nUpload aborted.");
		}
		extract_sequence_from_genome(\@lines,$filename,$genome);
	}
	my %hash=();
	open(OUT, ">$PATH_INFO/$loginname-config1.txt");
	open(INFO, "$PATH_INFO/$loginname-config.txt");
	while(my $line = <INFO>){
		if($file_type =~ /coordinates|attributes|conservation/){
			if($line =~ /^type_fasta=$fasta_file\s/){
				read_config_line($line,\%hash);
			}else{
				print OUT $line;
				}
		}elsif($file_type =~ /subset/){
			if($line =~ /^type_fasta=$fasta_file\s/){
				read_config_line($line,\%hash);
			}
			if($line !~ /^type_fasta=$filename\s/){
				print OUT $line;
			}
		}else{
			if($line =~ /^type_$file_type=$filename\s/){
				read_config_line($line,\%hash);
			}else{
				print OUT $line;
			}
		}
	}
	close INFO;
	if($file_type =~ /coordinates|attributes|conservation/){
		if(!%hash){ error_message("Fasta file not found!"); }
		$hash{$file_type} = $fasta_file;
		if($file_type eq "coordinates"){ $hash{$file_type} =~ s/\.fa$/\.coord/; }
		elsif($file_type eq "attributes"){ $hash{$file_type} =~ s/\.fa$/\.attr/; }
		elsif($file_type eq "conservation"){ $hash{$file_type} =~ s/\.fa$/\.cons/; }
		print OUT "type_fasta=$fasta_file";
		foreach my $key (keys %hash){
			if(!$key || $key eq "type_fasta"){ next; }
			print OUT "\t$key=$hash{$key}";
		}
	}elsif($file_type =~ /subset/){
		my $fileFasta = $fasta_file;
		if($fasta_file !~ /^public-/){
			$fileFasta = "$loginname-$fasta_file";
		}else{
			open(INFO, "$PATH_INFO/public-config.txt");
			my $fileFasta1 = $fasta_file;
			$fileFasta1 =~ s/^public-//;
			while(my $line = <INFO>){
				if($line =~ /^type_fasta=$fileFasta1/){
					read_config_line($line,\%hash); last;
				}
			}
		}
		$fileFasta = "$PATH_DATA/$fileFasta";
		if(subset_fasta($fileFasta,"$PATH_DATA/$loginname-$filename",\@lines)){
			my $description = $hashInput{"description"};
			print OUT "type_fasta=$filename\tdescription=$description";
			foreach my $key (keys %hash){
				my $newFilename = $filename;
				my $sourceFile = $fileFasta;
				if($key eq "coordinates"){ $newFilename =~ s/\.fa$/\.coord/; $sourceFile =~ s/\.fa$/\.coord/; }
				elsif($key eq "attributes"){ $newFilename =~ s/\.fa$/\.attr/; $sourceFile =~ s/\.fa$/\.attr/; }
				elsif($key eq "conservation"){ $newFilename =~ s/\.fa$/\.cons/; $sourceFile =~ s/\.fa$/\.cons/; }
				else{ next; }
				print OUT "\t$key=$newFilename";
				if($key =~ /coordinates|attributes/){
					subset_table($sourceFile,"$PATH_DATA/$loginname-$newFilename",\@lines);
				}elsif($key =~ /conservation/){
					subset_fasta($sourceFile,"$PATH_DATA/$loginname-$newFilename",\@lines);
				}
			}
		}
	}elsif($file_type =~ /subregion/){
		file_split($fasta_file,$filename,$description,\@lines);
	}else{
		print OUT "type_$file_type=$filename";
		if($description){ $hash{"description"} = $description; }
		foreach my $key (keys %hash){
			if(!$key || $key eq "type_$file_type"){ next; }
			print OUT "\t$key=$hash{$key}";
		}
	}
	print OUT "\n";
	close OUT;
	system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
	system("rm $PATH_INFO/$loginname-config1.txt");
	if($file_type !~ /subset|subregion/){
		if($file_type =~ /coordinates|attributes|conservation/){
			$filename = $hash{$file_type};
		}
		open(OUT, ">$PATH_DATA/$loginname-$filename");
		print OUT $list;
		close OUT;
	}
	if($file_type =~ /fasta|subset|subregion/){ $hashInput{"file_test"} = $filename; }
	elsif($file_type =~ /motif|pattern/){ $hashInput{"file_motif1"} = $filename; }
	elsif($file_type eq "search"){ $hashInput{"file_search"} = $filename; }
	elsif($file_type eq "repeat"){ $hashInput{"file_repeat"} = $filename; }
}else{
	print "content-type: text/html","\n\n";
	my @items = split(/[\&=]/,$form_info);
	while(@items){
		my $key = shift(@items);
		my $value = shift(@items);
		$key = substitute_chars($key);
		$value = substitute_chars($value);
		$value =~ s/^\s*//; $value =~ s/\s*$//;
		$key =~ s/^\s*//; $key =~ s/\s*$//;
		if($key eq ""){
			error_message("Empty key!","continue");
		}
		if($key=~/</ || $value=~/</){ error_message("Hyperlinks are blocked!"); }
		$hashInput{$key} = $value;
	}
	$loginname = $hashInput{"loginname"};
	$passwd = $hashInput{"passwd"};
	if($hashInput{"register"}){
		register_new_user();
	}elsif($hashInput{"login"}){
		if($loginname =~ /^guest$/i){
			my $xxx = int(rand 1000);
			$loginname = "guest$date$xxx";
			$hashInput{"loginname"} = $loginname;
			`cp $PATH_INFO/guest-config.txt $PATH_INFO/$loginname-config.txt`;
		}
		validate($loginname,$passwd);
		$sessionID = start_session($loginname);
	}elsif($hashInput{"reset_password"} && $loginname){
		$sessionID = $hashInput{"sessionID"};
		check_sessionID($loginname,$sessionID);
		reset_password($loginname);
	}else{
		$sessionID = $hashInput{"sessionID"};
		check_sessionID($loginname,$sessionID);
	}
}
#exit(0);
#print "FORM: $form_info\n";
clean_up();
my $dataname = $hashInput{"data"};
my $action = $hashInput{"action"};
if($action eq "change_password"){
	change_password_form($loginname,$sessionID);
}
elsif($action eq "update_password"){
	$passwd = $hashInput{"passwd"};
	validate($loginname,$passwd);
	update_password($loginname);
	$action = "continue";
}
my $file_delete = $hashInput{"file_delete"};
my $file_motif = $hashInput{"file_motif1"};
my $file_test = $hashInput{"file_test"};
#print "$dataname - $action - $file_motif - $file_test<br>\n";
if(!file_exist("$PATH_INFO/$loginname-config.txt")){
	`cp $PATH_INFO/default-config.txt $PATH_INFO/$loginname-config.txt`;
}
if($hashInput{"mainPage"} && $action =~ /continue|motif_get|motif_search|motif_improve|motif_compare|motif_cluster/){
	check_configuration();
}
######################### MAIN SUBROUTINE SWITCH ####################################
#if($hashInput{"update"}){				project_update(); }
if($hashInput{"motif_delete"}){				motif_delete(); }
if($hashInput{"motif_reverse"}){			motif_reverse(); }
if($hashInput{"file_delete"}){				file_delete(); }
if($hashInput{"file_public"}){				file_public(); }
if($hashInput{"action"} eq "save_motifs"){		save_motifs(); }
if($hashInput{"action"} eq "save_search"){		save_search(); }
#if($hashInput{"action"} eq "save_page"){		save_page(); }
#if($hashInput{"action"} eq "project_delete"){		project_delete(); $hashInput{"action"}="continue"; }
if($hashInput{"toregister"}){				print_registration_form(); }
elsif($hashInput{"file_download"}){			file_download(); }
elsif($hashInput{"action"} eq "generate_PFM"){		generate_PFM(); }
elsif($hashInput{"action"} eq "motif_get"){		motif_get(); }
elsif($hashInput{"action"} eq "motif_get1"){		motif_get1(); }
elsif($hashInput{"action"} eq "motif_search"){		motif_search(); }
elsif($hashInput{"action"} eq "motif_search1"){		motif_search1(); }
elsif($hashInput{"action"} eq "motif_improve"){		motif_improve(); }
elsif($hashInput{"action"} eq "motif_improve1"){	motif_improve1(); }
elsif($hashInput{"action"} eq "motif_compare"){		motif_compare(); }
elsif($hashInput{"action"} eq "motif_compare1"){	motif_compare1(); }
elsif($hashInput{"action"} eq "motif_cluster"){		motif_cluster(); }
elsif($hashInput{"action"} eq "motif_cluster1"){	motif_cluster1(); }
#elsif($hashInput{"action"} eq "project_new"){		project_edit(); }
#elsif($hashInput{"action"} eq "project_edit"){		project_edit($dataname); }
elsif($hashInput{"action"} eq "save_PFM"){		save_PFM(); }
elsif($hashInput{"action"} eq "search_motif_name"){	search_motif_name(); }
elsif($hashInput{"action"} eq "show_cluster"){		show_cluster(); }
elsif($hashInput{"action"} eq "show_motifs"){		show_motifs(); }
elsif($hashInput{"action"} eq "show_motif"){		show_motif(); }
#elsif($hashInput{"action"} eq "show_project"){		show_project(); }
elsif($hashInput{"action"} eq "show_search_results"){	show_search_results(); }
elsif($hashInput{"action"} eq "show_search_results1"){	show_search_results1(); }
elsif($hashInput{"action"} eq "show_search_motifs"){	show_search_motifs(); }
elsif($hashInput{"action"} eq "show_search_motif"){	show_search_motif(); }
elsif($hashInput{"action"} eq "show_search_sequences"){	show_search_sequences(); }
elsif($hashInput{"action"} eq "show_search_sequence"){	show_search_sequence(); }
elsif($hashInput{"action"} eq "show_search_sequence_seq"){show_search_sequence_seq(); }
elsif($hashInput{"login"} || $hashInput{"action"} eq "continue"){ main_page(); }
else{ 							error_message("Unknown command!"); }
exit(0);

#************************************
sub  main_page
#************************************
{
my @data_list;
my @saved_pages;
my @motif_list;
my @fasta_list;
my @repeat_list;
my @search_list;
my %hashDefault;
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
while(my $line = <INFO>){
	my %hash=();
	if(length($line) < 3) { next; }
	read_config_line($line,\%hash);
	my $saved_page = $hash{"saved_page"};
	if($saved_page){ push(@saved_pages,[$saved_page,$hash{"description"}]); next; }
	my $file_motif = $hash{"type_motif"};
	if($file_motif){ push(@motif_list,[$file_motif,$hash{"description"}]); next; }
	my $file_repeat = $hash{"type_repeat"};
	if($file_repeat){ push(@repeat_list,[$file_repeat,$hash{"description"}]); next; }
	my $file_fasta = $hash{"type_fasta"};
	if($file_fasta){ push(@fasta_list,[$file_fasta,$hash{"description"}]); next; }
	my $file_search = $hash{"type_search"};
	if($file_search){ push(@search_list,[$file_search,$hash{"description"}]); next; }
	my $data = $hash{"data"};
	if($hash{"mainPage"}){ %hashDefault = %hash; next; }
	elsif($data){
		push(@data_list, [$data,$hash{"description"},$hash{"file_test"}]);
	}
}
close INFO;
@saved_pages = sort {lc($a->[0]) cmp lc($b->[0])} @saved_pages;
@motif_list = sort {lc($a->[0]) cmp lc($b->[0])} @motif_list;
@repeat_list = sort {lc($a->[0]) cmp lc($b->[0])} @repeat_list;
@fasta_list = sort {lc($a->[0]) cmp lc($b->[0])} @fasta_list;
@search_list = sort {lc($a->[0]) cmp lc($b->[0])} @search_list;
@data_list = sort {lc($a->[0]) cmp lc($b->[0])} @data_list;
if($loginname ne "public" && open(INFO, "$PATH_INFO/public-config.txt")){
	my @motif_list1=();
	my @fasta_list1=();
	my @search_list1=();
	my @repeat_list1=();
	while(my $line = <INFO>){
		my %hash=();
		read_config_line($line,\%hash);
		my $file_motif = $hash{"type_motif"};
		if($file_motif){ push(@motif_list1,["public-$file_motif",$hash{"description"}]); next; }
		my $file_repeat = $hash{"type_repeat"};
		if($file_repeat){ push(@repeat_list1,["public-$file_repeat",$hash{"description"}]); next; }
		my $file_fasta = $hash{"type_fasta"};
		if($file_fasta){ push(@fasta_list1,["public-$file_fasta",$hash{"description"}]); }
		my $file_search = $hash{"type_search"};
		if($file_search){ push(@search_list1,["public-$file_search",$hash{"description"}]); }
	}
	push(@motif_list, sort {lc($a->[0]) cmp lc($b->[0])} @motif_list1);
	push(@repeat_list, sort {lc($a->[0]) cmp lc($b->[0])} @repeat_list1);
	push(@fasta_list, sort {lc($a->[0]) cmp lc($b->[0])} @fasta_list1);
	push(@search_list, sort {lc($a->[0]) cmp lc($b->[0])} @search_list1);
}
close INFO;
my $interval = $hashDefault{"interval"};
if(!$interval){ $interval=100; }
my $match_thresh = $hashDefault{"match_thresh"};
if(!$match_thresh){ $match_thresh=0.75; }

print "<HTML><HEAD><TITLE>CisFinder</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
my ($items,$descriptions,$files) = get_array_lists(\@data_list);
print "data_list = new Array($items);\n";
print "data_description = new Array($descriptions);\n";
print "data_files = new Array($files);\n";
($items,$descriptions) = get_array_lists(\@motif_list);
print "motif_list = new Array($items);\n";
print "motif_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@repeat_list);
print "repeat_list = new Array($items);\n";
print "repeat_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@fasta_list);
print "fasta_list = new Array($items);\n";
print "fasta_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@search_list);
print "search_list = new Array($items);\n";
print "search_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@saved_pages);
print "saved_page_list = new Array($items);\n";
print "saved_page_description = new Array($descriptions);\n";
print "function update_description() {\n";
print "  var index;\n";
if($loginname !~ /^guest/){
	if(@data_list){
		print "  index = document.cisfinder.data.selectedIndex;\n";
		print "  document.cisfinder.description_data.value = data_description[index];\n";
	}
	if(@saved_pages){
		print "  index = document.cisfinder.saved_page.selectedIndex;\n";
		print "  document.cisfinder.description_saved_page.value = saved_page_description[index];\n";
	}
}
print "  index = document.cisfinder.file_test.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_test.value = fasta_description[index-1]; }\n";
print "  else{ document.cisfinder.description_test.value = \"\"; }\n";
print "  index = document.cisfinder.file_control.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_control.value = fasta_description[index-1]; }\n";
print "  else{ document.cisfinder.description_control.value = \"\"; }\n";
print "  index = document.cisfinder.file_repeat.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_repeat.value = repeat_description[index-1]; }\n";
print "  else{ document.cisfinder.description_repeat.value = \"\"; }\n";
print "  index = document.cisfinder.file_motif1.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_motif1.value = motif_description[index-1]; }\n";
print "  else{ document.cisfinder.description_motif1.value = \"\"; }\n";
print "  index = document.cisfinder.file_search.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_search.value = search_description[index-1]; }\n";
print "  else{ document.cisfinder.description_search.value = \"\"; }\n";
print "}\n";
print "function clear_marks() {\n";
print "	document.cisfinder.target = \"\"\n";
print "	document.cisfinder.file_download.value = \"\";\n";
print "	document.cisfinder.file_delete.value = \"\";\n";
print "	document.cisfinder.file_public.value = \"\";\n";
print "	document.cisfinder.mainPage.value = \"mainPage\";\n";
print "}\n";
print "function do_analysis(command) {\n";
print "	if(command==\"motif_search\" || command==\"motif_improve\" || command==\"motif_cluster\" || command==\"motif_compare\" || command==\"show_motifs\"){\n";
print "		if(document.cisfinder.file_motif1.selectedIndex==0){\n";
print "			alert(\"Select Motif file!\"); return false;\n";
print "		}\n";
print "	}\n";
print "	if(command==\"show_search_results\"){\n";
print "		if(document.cisfinder.file_search.selectedIndex==0){ alert(\"File with search results not selected\"); return false; }\n";
print "	}\n";
print "	clear_marks();\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	if(command == \"show_project\"){\n";
print "		document.cisfinder.mainPage.value = \"\";\n";
print "	}\n";
print "	else{\n";
print "		document.cisfinder.target = \"_BLANK\"+x;\n";
print "	}\n";
print "	document.cisfinder.action.value = command;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
#print "function data_manage(command){\n";
#print "	if(command==\"project_delete\" && data_list.length<=1){ alert(\"Last data set cannot be deleted\"); return false; }\n";
#if($loginname =~ /^guest/){
#	print "	if(command==\"project_delete\" || command==\"project_new\"){\n";
#	print "		alert(\"Command not available for guests\"); return false;\n";
#	print "	}\n";
#}
#if($loginname !~ /^public|^administrator/){
#	print "var data = document.cisfinder.data.options[document.cisfinder.data.selectedIndex].value;\n";
#	print "if(data.search(/^public-/) >=0){ alert(\"Public projects cannot be deleted\"); return false; }\n";
#}
#print "	clear_marks();\n";
#print "	document.cisfinder.action.value = command;\n";
#print "	document.cisfinder.submit();\n";
#print "}\n";
print "function download_file(file_type) {\n";
print "	var filename;\n";
print "	if(file_type == \"file_fasta\"){ filename = document.cisfinder.select_fasta_file.options[document.cisfinder.select_fasta_file.selectedIndex].value; }\n";
print "	else if(file_type == \"file_repeat\"){ filename = document.cisfinder.select_repeat_file.options[document.cisfinder.select_repeat_file.selectedIndex].value; }\n";
print "	else if(file_type == \"file_motif\"){ filename = document.cisfinder.select_motif_file.options[document.cisfinder.select_motif_file.selectedIndex].value; }\n";
print "	else if(file_type == \"file_search\"){ filename = document.cisfinder.select_search_file.options[document.cisfinder.select_search_file.selectedIndex].value; }\n";
print "	if(!filename){ alert(\"File not selected\"); return false; }\n";
print "	clear_marks();\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"\";\n";
print "	document.cisfinder.file_download.value = filename;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function edit_refresh() {\n";
print "	clear_marks();\n";
print "	document.cisfinder.action.value = \"continue\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function upload_onsubmit() {\n";
print "	document.cisfinder.action.value = \"continue\";\n";
print "	if(document.upload.upload_file_name.value && document.upload.pasted_text.value){\n";
print "		alert(\"You cannot use both text and file for submission\\nRemove one of them\"); return false;\n";
print "	}\n";
print "	if(!document.upload.upload_file_name.value && !document.upload.pasted_text.value){\n";
print "		if(document.upload.file_type.value==\"subregion\"){\n";
print "			alert(\"Enter coordinate region(s) to extract\\n(e.g. '400-600' or '0-500;1500-2000')\");\n";
print "		}else{\n";
print "			alert(\"Paste text or select file to upload\");\n";
print "		}\n";
print "		return false;\n";
print "	}\n";
print "	var file = document.upload.upload_file_name.value;\n";
print "	if(document.upload.rename.value){ file=document.upload.rename.value; }\n";
print "	else if(document.upload.pasted_text.value){\n";
print "		if(document.upload.file_type.value==\"fasta\"){ file=\"user_file.fa\"; }\n";
print "		else if(document.upload.file_type.value==\"motif\" || document.upload.file_type.value==\"pattern\"){ file=\"user_motif\"; }\n";
print "		else if(document.upload.file_type.value==\"search\"){ file=\"user_search\"; }\n";
print "		else if(document.upload.file_type.value==\"repeat\"){ file=\"user_repeat\"; }\n";
print "		else{ file=\"user_file\"; }\n";
print "	}\n";
print "	if(document.upload.file_type.value==\"fasta\" && file.search(/\\.fa\$/) < 0){\n";
print "		alert(\"NOTE: Fasta file name will be appended '.fa' extension\");\n";
print "	}\n";
print "	var i;\n";
print "	for(i=0; i<fasta_list.length; ++i){\n";
print "		if(file == fasta_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	for(i=0; i<repeat_list.length; ++i){\n";
print "		if(file == repeat_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	for(i=0; i<motif_list.length; ++i){\n";
print "		if(file == motif_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	for(i=0; i<search_list.length; ++i){\n";
print "		if(file == search_list[i] && !confirm(\"File \"+file+\" already exists. Do you want to overwrite it?\")){\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	if(file.search(/ |\\/|\\&|\\t|\\%|\\+|\\@|\\!/) >= 0){\n";
print "		alert(\"File name should have neither spaces nor special characters.\\nRename it provided field\");\n";
print "		return(false);\n";
print "	}\n";
print "	if(document.upload.file_type.selectedIndex>4){\n";
print "		if(document.upload.file_type.selectedIndex==5){\n";
print "			if(!document.upload.fasta_file.selectedIndex && !confirm(\"You have not selected the associated fasta file.\\nDo you want to extract sequence from genome (mm9, hg19, or hg18)?\")){\n";
print "				return false;\n";
print "			}\n";
print "		}\n";
print "		else if(!document.upload.fasta_file.selectedIndex){\n";
print "			alert(\"Select the associated fasta file\");\n";
print "			return(false);\n";
print "		}\n";
print "		if(document.upload.file_type.selectedIndex>=8 || (document.upload.file_type.selectedIndex==5 && !document.upload.fasta_file.selectedIndex)){\n";
print "			file = document.upload.rename.value;\n";
print "			if(!file){\n";
print "				alert(\"Enter new file name (fasta) with '.fa' ending.\");\n";
print "				return false;\n";
print "			}\n";
print "			if(file.search(/\\.fa\$/) < 0){\n";
print "				alert(\"NOTE: Fasta file name will be appended '.fa' extension\");\n";
print "			}\n";
print "		}\n";
print "	}\n";
print "}\n";
print "function delete_file(file_type) {\n";
print "	var filename;\n";
print "	if(file_type == \"fasta\"){\n";
print "		filename = document.cisfinder.select_fasta_file.options[document.cisfinder.select_fasta_file.selectedIndex].value;\n";
print "		var i;\n";
print "		for(i=0; i<data_files.length; ++i){\n";
print "			if(filename == data_files[i]){\n";
print "				alert(\"File is used in project \"+data_list[i]+\". First remove the file from project.\"); return false;\n";
print "			}\n";
print "		}\n";
print "		document.cisfinder.select_fasta_file.selectedIndex = 0;\n";
print "	}\n";
print "	else if(file_type == \"repeat\"){ filename = document.cisfinder.select_repeat_file.options[document.cisfinder.select_repeat_file.selectedIndex].value; }\n";
print "	else if(file_type == \"motif\"){ filename = document.cisfinder.select_motif_file.options[document.cisfinder.select_motif_file.selectedIndex].value; }\n";
print "	else if(file_type == \"search\"){ filename = document.cisfinder.select_search_file.options[document.cisfinder.select_search_file.selectedIndex].value; }\n";
print "	if(!filename){ alert(\"File not selected\"); return false; }\n";
if($loginname !~ /^public/){
	print "if(filename.search(/^public-/) >=0){ alert(\"Public files cannot be deleted\"); return false; }\n";
}
print "	clear_marks();\n";
print "	document.cisfinder.file_delete.value = filename;\n";
print "	document.cisfinder.action.value = \"continue\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function explore_motifs() {\n";
print "	if(!document.cisfinder.file_motif1.options[document.cisfinder.file_motif1.selectedIndex].value){ alert(\"Motif file not selected\"); return false; }\n";
print "	clear_marks();\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.action.value = \"show_motifs\";\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function open_page() {\n";
print "	var source = \"$HOME_ADDRESS/saved/$loginname-\" + document.cisfinder.saved_page.options[document.cisfinder.saved_page.selectedIndex].value + \".html\";\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	var name = \"_BLANK\"+x;\n";
print "	window.open(source,name);\n";
print "}\n";
print "function delete_page() {\n";
if($loginname =~ /^guest/){ print "alert(\"Command not available for guests\"); return false;\n"; }
print "	document.cisfinder.file_delete.value = \"saved/$loginname-\" + document.cisfinder.saved_page.options[document.cisfinder.saved_page.selectedIndex].value;\n";
print "	document.cisfinder.action.value = \"continue\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function select_string(string,select_object) {\n";
print "	var i;\n";
print "	for(i=0;i<select_object.options.length && string != select_object.options[i].value; ++i){}\n";
print "	if(i==select_object.options.length){ i=0; }\n";
print "	select_object.selectedIndex = i;\n";
print "}\n";
print "function run_demo() {\n";
print "	select_string(\"public-Pou5f1_binding.fa\",document.cisfinder.file_test);\n";
print "	select_string(\"public-Pou5f1_control.fa\",document.cisfinder.file_control);\n";
print "	select_string(\"public-motif_Pou5f1\",document.cisfinder.file_motif1);\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header("update_description();");
print "HELP: <a href=../cisfinder-help.html TARGET=\"_blank324\">A guide to CisFinder</a>\n";
print "(if you don't see help <a href=../cisfinder-help.html#scripted>enable scripted windows</a>)<br>\n";
print "<FORM NAME=upload ENCTYPE=\"multipart/form-data\" ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST onSubmit=\"return upload_onsubmit();\">\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"mainPage\" TYPE=\"hidden\" VALUE=\"mainPage\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"continue\">\n";
print "<FONT SIZE=+1><b>Upload New Data</b></FONT> Paste or use file; check <a href=../cisfinder-help.html#format>file format</a> before uploading!<br>";
print "<TEXTAREA NAME=pasted_text ROWS=3 COLS=80></TEXTAREA><br>\n";
print "Use file:<INPUT NAME=\"upload_file_name\" TYPE=\"file\">\n";
print "File type: <SELECT NAME=file_type>";
print "<OPTION VALUE=fasta>Fasta";
print "<OPTION VALUE=motif>Motif (PFM)";
print "<OPTION VALUE=pattern>Motif (consensus)";
print "<OPTION VALUE=repeat>Repeat";
print "<OPTION VALUE=search>Search results";
print "<OPTION VALUE=coordinates>Genome location";
print "<OPTION VALUE=attributes>Sequence attributes";
print "<OPTION VALUE=conservation>Conservation\n";
print "<OPTION VALUE=subset>Subset of sequences\n";
print "<OPTION VALUE=subregion>Sub-regions of sequences\n";
print "</SELECT> &nbsp; \n";
print "<input type=\"submit\" name=\"load_file\" value=\"Upload\"><br>\n";
print "<TABLE BORDER=0><TR><TD>Rename file as:<TD><INPUT NAME=\"rename\" SIZE=30> (<i>Optional</i>) No spaces; use '.fa' for fasta\n";
print "<TR><TD>Description:<TD><INPUT NAME=\"description\" SIZE=80>\n";
print "<TR><TD>Associate with:\n";
print "<TD><select name=\"fasta_file\">\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@fasta_list; ++$i){ 
	print "<option value=\"$fasta_list[$i]->[0]\"> $fasta_list[$i]->[0]\n";
}
print "</select> (For special file types)\n";
print "</TABLE>\n";
print "</FORM>\n";
print "<HR NOSHADE></HR>\n";
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST LANGUAGE=\"javascript\" onsubmit=\"return cisfinder_onsubmit()\">\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"continue\">\n";
print "<INPUT NAME=\"file_download\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"file_delete\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"file_public\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"mainPage\" TYPE=\"hidden\" VALUE=\"mainPage\">\n";
print "<INPUT NAME=\"interval\" TYPE=\"hidden\" VALUE=\"$interval\">\n";
print "<TABLE border=0>\n";

print "<TR><TD><INPUT NAME=\"button_refresh\" TYPE=\"button\" VALUE=\"Refresh\" onClick=\"edit_refresh();\"> \n";
print "<FONT SIZE=+1><b>Select Data Files<TD><b><center>Description<TD><b>File type\n";
print "<TR><TD><select name=\"file_test\" onChange=update_description();>\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@fasta_list; ++$i){ 
	print "<option value=\"$fasta_list[$i]->[0]\"";
	if($fasta_list[$i]->[0] eq $hashDefault{"file_test"}){ print " selected"; }
	print "> $fasta_list[$i]->[0]\n";
}
print "</select><td><INPUT NAME=\"description_test\" SIZE=40><TD>Sequence file #1 (test)\n";
print "<TR><TD><select name=\"file_control\" onChange=update_description();>\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@fasta_list; ++$i){ 
	print "<option value=\"$fasta_list[$i]->[0]\"";
	if($fasta_list[$i]->[0] eq $hashDefault{"file_control"}){ print " selected"; }
	print "> $fasta_list[$i]->[0]\n";
}
print "</select><td><INPUT NAME=\"description_control\" SIZE=40><TD> Sequence file #2 (control)\n";
print "<TR><TD><select name=\"file_motif1\" onChange=update_description();>\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@motif_list; ++$i){ 
	print "<option value=\"$motif_list[$i]->[0]\"";
	if($motif_list[$i]->[0] eq $hashDefault{"file_motif1"}){ print " selected"; }
	print "> $motif_list[$i]->[0]\n";
}
print "</select><td><INPUT NAME=\"description_motif1\" SIZE=40><TD>Motif file\n";
print "<TR><TD><select name=\"file_repeat\" onChange=update_description();>\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@repeat_list; ++$i){ 
	print "<option value=\"$repeat_list[$i]->[0]\"";
	if($repeat_list[$i]->[0] eq $hashDefault{"file_repeat"}){ print " selected"; }
	print "> $repeat_list[$i]->[0]\n";
}
print "</select><td><INPUT NAME=\"description_repeat\" SIZE=40><TD>Repeat file\n";
print "<TR><TD><select name=\"file_search\" onChange=update_description();>\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@search_list; ++$i){ 
	print "<option value=\"$search_list[$i]->[0]\"";
	if($search_list[$i]->[0] eq $hashDefault{"file_search"}){ print " selected"; }
	print "> $search_list[$i]->[0]\n";
}
print "</select><td><INPUT NAME=\"description_search\" SIZE=40><TD>Search results\n";
print "</TABLE>\n";
print "<INPUT NAME=\"button_demo\" TYPE=\"button\" VALUE=\"  Demo  \" onClick=\"run_demo();\"> \n";
print "Click \"Demo\" and then \"Identify motifs\" (below). <FONT SIZE=-1>(Also, click on pull-down menus to see other public data sets available for analysis)</FONT><p>\n";

print "<FONT SIZE=+1><b>Analysis Tools</b></FONT> (Parameters are specified later)<br>\n";
print "<TABLE border=0>\n";
print "<TR><TD><INPUT NAME=\"get_motif_button\" TYPE=\"button\"      VALUE=\"Identify motifs\" style=width:200px; onClick=\"do_analysis('motif_get');\">\n";
print "<TD>Identify motifs over-repretented in sequence file #1 compared to sequence file #2 (#2 is optional)\n";
print "<TR><TD><INPUT NAME=\"improve_motif_button\" TYPE=\"button\"  VALUE=\"Improve motifs\" style=width:200px; onClick=\"do_analysis('motif_improve');\">\n";
print "<TD>Improve motifs from motif file based on occurrence in sequence files #1 and #2 (#2 is optional)\n";
print "<TR><TD><INPUT NAME=\"cluster_motif_button\" TYPE=\"button\"  VALUE=\"Cluster motifs\" style=width:200px; onClick=\"do_analysis('motif_cluster');\">\n";
print "<TD>Cluster motifs from motif file based on similarity\n";
print "<TR><TD><INPUT NAME=\"compare_motif_button\" TYPE=\"button\"  VALUE=\"Compare motifs\" style=width:200px; onClick=\"do_analysis('motif_compare');\">\n";
print "<TD>Compare motifs from motif file with other motif files (e.g. for annotation)\n";
print "<TR><TD><INPUT NAME=\"explore_motif_button\" TYPE=\"button\"  VALUE=\"Show motifs\" style=width:200px; onClick=\"explore_motifs();\">\n";
print "<TD>Shows motif logo and PFM for motif file (includes motif copying)\n";
print "<TR><TD><INPUT NAME=\"search_motif_button\" TYPE=\"button\"   VALUE=\"Search motifs\" style=width:200px; onClick=\"do_analysis('motif_search');\">\n";
print "<TD>Search sequence file #1 for occurence of motifs from motif file\n";
print "<TR><TD><INPUT NAME=\"search_results_button\" TYPE=\"button\" VALUE=\"Show search results\" style=width:200px; onClick=\"do_analysis('show_search_results');\">\n";
print "<TD>Show search results (motifs, sequences, frequency distribution)\n";
print "</TABLE><p>\n";

print "<FONT SIZE=+1><b>Data Management</b></FONT><br>\n";
print "<TABLE border=0>\n";
print "<TR><TD><select name=\"select_fasta_file\">\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@fasta_list; ++$i){ 
	print "<option value=\"$fasta_list[$i]->[0]\"> $fasta_list[$i]->[0]\n";
}
print "</select><TD>Sequence file\n";
print "<td><INPUT NAME=\"file_fasta_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_fasta');\">\n";
print "<td><INPUT NAME=\"delete_fasta_button\" TYPE=\"button\" VALUE=\"Delete\" onClick=\"delete_file('fasta');\">\n";
print "<TR><TD><select name=\"select_motif_file\">\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@motif_list; ++$i){ 
	print "<option value=\"$motif_list[$i]->[0]\"> $motif_list[$i]->[0]\n";
}
print "</select><TD>Motif file\n";
print "<td><INPUT NAME=\"file_motif_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_motif');\">\n";
print "<td><INPUT NAME=\"delete_motif_button\" TYPE=\"button\" VALUE=\"Delete\" onClick=\"delete_file('motif');\">\n";
print "<TR><TD><select name=\"select_repeat_file\">\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@repeat_list; ++$i){ 
	print "<option value=\"$repeat_list[$i]->[0]\"> $repeat_list[$i]->[0]\n";
}
print "</select><TD>Repeat file\n";
print "<td><INPUT NAME=\"file_repeat_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_repeat');\">\n";
print "<td><INPUT NAME=\"delete_repeat_button\" TYPE=\"button\" VALUE=\"Delete\" onClick=\"delete_file('repeat');\">\n";
print "<TR><TD><select name=\"select_search_file\">\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@search_list; ++$i){ 
	print "<option value=\"$search_list[$i]->[0]\"> $search_list[$i]->[0]\n";
}
print "</select><TD>Search file\n";
print "<td><INPUT NAME=\"file_search_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_search');\">\n";
print "<td><INPUT NAME=\"delete_search_button\" TYPE=\"button\" VALUE=\"Delete\" onClick=\"delete_file('search');\">\n";
print "</TABLE><p>\n";
#if($loginname !~ /^guest/){
#	print "<TABLE border=0>\n";
#	if(@saved_pages){
#		print "<tr><td><FONT SIZE=+1><b>Saved page<td><b><center>Description\n";
#		print "<tr><td><SELECT NAME=\"saved_page\" onChange=update_description();>\n";
#		foreach my $ref (@saved_pages){
#			my $page_name = $ref->[0];
#			print "<OPTION VALUE=$page_name> $page_name\n";
#		}
#		print "</SELECT>\n";
#		print "<td><INPUT NAME=\"description_saved_page\" SIZE=58>\n";
#		print "<td><INPUT TYPE=button name=open_page_button value=\" Open page  \" onClick=open_page();>\n";
#		if($loginname !~ /^guest/){
#			print "<td><INPUT TYPE=button name=delete_page_button value=\"Delete\"  onClick=delete_page();>\n";
#		}
#	}
#	if(@data_list){
#		print "<tr><td><FONT SIZE=+1><b>Project<td><b><center>Description\n";
#		print "<tr><td><select name=\"data\" onChange=update_description();>\n";
#		for(my $i=0; $i<@data_list; ++$i){ 
#			print "<option value=\"$data_list[$i]->[0]\">$data_list[$i]->[0]";
#		}
#		print "</select>\n";
#		print "<td><INPUT NAME=\"description_data\" SIZE=58>\n";
#		print "<td><INPUT NAME=\"open_data_button\" TYPE=\"button\" VALUE=\"Open project\" onClick=\"do_analysis('show_project');\">\n";
#		print "</TABLE>\n";
#	}
#	print "<INPUT NAME=\"new_data_button\" TYPE=\"button\" VALUE=\"New project\" onClick=\"data_manage('project_new');\">\n";
#	if(@data_list){
#		print "<INPUT NAME=\"project_delete\" TYPE=\"button\" VALUE=\"Delete project\" onClick=\"data_manage('project_delete');\">\n";
#	}
#}
print "</FORM>\n";
my $x = int(10000*rand());
print "<p><FORM NAME=change_pas ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST TARGET=\"_blank$x\">\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=change_password>\n";
print "<INPUT TYPE=\"submit\" VALUE=\"   Change password   \">\n";
print "</table></FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#************************************
sub  show_project
#************************************
{
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
my %hashProject;
my %fileDescription;
my @saved_pages=();
while(my $line = <INFO>){
	my %hash=();
	if(length($line) < 3) { next; }
	read_config_line($line,\%hash);
	if($hash{"mainPage"}){ next; }
	my $data = $hash{"data"};
	if($data eq $dataname){ %hashProject = %hash; }
	my $page = $hash{"saved_page"};
	if($page){ push(@saved_pages,$page); }
}
close INFO;
my $file_test = $hashProject{"file_test"};
my $file_control = $hashProject{"file_control"};
my $file_repeat = $hashProject{"file_repeat"};
my $file_motif1 = $hashProject{"file_motif1"};
my $file_search = $hashProject{"file_search"};
my $description = $hashProject{"description"};
my $use_repeats = $hashProject{"use_repeats"};
my $FDR = $hashProject{"FDR"};
my $strand = $hashProject{"strand"};
my $score = $hashProject{"score"};
my $adjust_cg = $hashProject{"adjust_cg"};
my $presence = $hashProject{"presence"};
my $min_ratio = $hashProject{"min_ratio"};
my $max_repeat = $hashProject{"max_repeat"};
my $match_thresh = $hashProject{"match_thresh"};
my $method = $hashProject{"method"};
my $interval = $hashProject{"interval"};
my $Nmatrix = $hashProject{"Nmatrix"};
if(!$Nmatrix){ $Nmatrix=500; }

my ($description_test,$description_control,$description_repeat,$description_motif1,$description_motif2,$description_search);
open(INFO, "$PATH_INFO/$loginname-config.txt");
while(my $line = <INFO>){
	my %hash=();
	if(length($line) < 3) { next; }
	read_config_line($line,\%hash);
	if($file_test && $file_test eq $hash{"type_fasta"}){ $description_test = $hash{"description"}; next; }
	if($file_control && $file_control eq $hash{"type_fasta"}){ $description_control = $hash{"description"}; next; }
	if($file_repeat && $file_repeat eq $hash{"type_repeat"}){ $description_repeat = $hash{"description"}; next; }
	if($file_motif1 && $file_motif1 eq $hash{"type_motif"}){ $description_motif1 = $hash{"description"}; next; }
	if($file_search && $file_search eq $hash{"type_search"}){ $description_search = $hash{"description"}; next; }
}
close INFO;

print "<HTML><HEAD><TITLE>CisFinder</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";

my $page_list = join("\",\"",@saved_pages);
if($page_list){ $page_list = "\"".$page_list."\""; }
print "saved_pages = new Array($page_list);\n";
print "function clear_marks() {\n";
print "	document.cisfinder.target = \"\"\n";
print "	document.cisfinder.file_download.value = \"\";\n";
print "}\n";
print "function do_analysis(command) {\n";
print "	clear_marks();\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = command;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function download_file(file_type) {\n";
print "	var filename;\n";
print "	if(file_type == \"file_test\"){ filename = document.cisfinder.file_test.value; }\n";
print "	if(file_type == \"file_control\"){ filename = document.cisfinder.file_control.value; }\n";
print "	if(file_type == \"file_motif1\"){ filename = document.cisfinder.file_motif1.value; }\n";
print "	if(file_type == \"file_repeat\"){ filename = document.cisfinder.file_repeat.value; }\n";
print "	if(file_type == \"file_search\"){ filename = document.cisfinder.file_search.value; }\n";
print "	if(!filename){ alert(\"File not selected\"); return false; }\n";
print "	clear_marks();\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"\";\n";
print "	document.cisfinder.file_download.value = filename;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function data_manage(command){\n";
print "	clear_marks();\n";
print "	document.cisfinder.action.value = command;\n";
print "	document.cisfinder.file_download.value = \"\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function explore_motifs() {\n";
print "	clear_marks();\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.action.value = \"show_motifs\";\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function cancel_project() {\n";
print "	clear_marks();\n";
print "	document.cisfinder.action.value = \"continue\";\n";
print "	document.cisfinder.data.value = \"\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
#print "function save_page(){\n";
#print "	var page_name = document.cisfinder.web_page_name.value;\n";
#print "	if(!page_name){ alert(\"You need to enter the name for web page\"); return false; }\n";
#print "	if(page_name.search(/ |\\/|\\&|\\t|\\%|\\+|\\@|\\!/) >= 0){\n";
#print "		alert(\"Web page name should have neither spaces nor special characters\");\n";
#print "		return(false);\n";
#print "	}\n";
#print "	var i;\n";
#print "	for(i=0; i<saved_pages.length; ++i){\n";
#print "		if(page_name == saved_pages[i]){\n";
#print "			if(!confirm(\"File already exists. Do you want to overwrite it?\")){ return false; }\n";
#print "		}\n";
#print "	}\n";
#print "	clear_marks();\n";
#print "	var x = Math.round(Math.random()*1000);\n";
#print "	document.cisfinder.target = \"_BLANK\"+x;\n";
#print "	document.cisfinder.action.value = \"save_page\";\n";
#print "	document.cisfinder.submit();\n";
#print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "HELP: <a href=\"$HOME_ADDRESS/cisfinder-help.html\" TARGET=\"_blank324\">A guide to CisFinder</a><br>\n";
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST LANGUAGE=\"javascript\" onsubmit=\"return cisfinder_onsubmit()\">\n";
#print "<H3>Project: <font color=brown>$dataname</font>\n";
#if($loginname !~ /^guest/){
#	print "<INPUT NAME=\"project_edit_button\" TYPE=\"button\" VALUE=\"Edit project\" onClick=\"data_manage('project_edit');\">\n";
#}
#print "<INPUT NAME=\"project_abort\" TYPE=\"button\" VALUE=\"Leave project\" onClick=\"cancel_project();\">\n";
print "</h3><INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"analysis\">\n";
print "<INPUT NAME=\"file_download\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<b>Description:</b> $description<p>\n";
print "<TABLE border=1>\n";
print "<tr><td><b>File type<td><b>File name<td><b>Description<td><b>Download\n";
print "<tr><td><b>Sequence file #1 (test)<td><font color=brown>$file_test<td>$description_test\n";
print "<INPUT NAME=\"file_test\" TYPE=\"hidden\" VALUE=\"$file_test\">\n";
print "<td><INPUT NAME=\"file_test_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_test');\">\n";
if($file_control){
	print "<tr><td><b>Sequence file #2 (control)<td><font color=brown>$file_control<td>$description_control\n";
	print "<td><INPUT NAME=\"file_control_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_control');\">\n";
	print "<INPUT NAME=\"file_control\" TYPE=\"hidden\" VALUE=\"$file_control\">\n";
}
if($file_repeat){
	print "<tr><td><b>Repeat file<td><font color=brown>$file_repeat<td>$description_repeat\n";
	print "<td><INPUT NAME=\"file_repeat_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_repeat');\">\n";
	print "<INPUT NAME=\"file_repeat\" TYPE=\"hidden\" VALUE=\"$file_repeat\">\n";
}
if($file_motif1){
	print "<tr><td><b>Motif file<td><font color=brown>$file_motif1<td>$description_motif1\n";
	print "<td><INPUT NAME=\"file_motif1_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_motif1');\">\n";
	print "<INPUT NAME=\"file_motif1\" TYPE=\"hidden\" VALUE=\"$file_motif1\">\n";
}
if($file_search){
	print "<tr><td><b>Search results<td><font color=brown>$file_search<td>$description_search\n";
	print "<td><INPUT NAME=\"file_search_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_search');\">\n";
	print "<INPUT NAME=\"file_search\" TYPE=\"hidden\" VALUE=\"$file_search\">\n";
}
print "</TABLE>\n";

print "<h3>Parameters</h3>\n";
my @yesno = ("No","Yes");
my @strand_list = ("Both","Forward","Forward vs. back");
my @score_list = ("Z","Z+ratio","Z+info","Z+ratio+info","Z+repeat","Z+ratio+repeat","Z+info+repeat","Z+ratio+info+repeat");
my @method_list = ("Similarity","Similarity and linkage","Linkage");
print "<table border=1>\n";
print "<tr><td><b>Parameter<td><b>Value<td><b>Parameter<td><b>Value\n";
print "<tr><td>Search repeats<td><font color=brown>$yesno[$use_repeats]\n";
print "<td>Search in strands<td><font color=brown>$strand_list[$strand]\n";
print "<tr><td>False Discovery Rate (FDR)<td><font color=brown>$FDR\n";
print "<td>Adjust for CG/AT ratio and CpG<td><font color=brown>$yesno[$adjust_cg]\n";
print "<tr><td>Count motif once per sequence<td><font color=brown>$yesno[$presence]\n";
print "<td>Score motifs by<td><font color=brown>$score_list[$score]\n";
print "<tr><td>Min enrichment ratio<td><font color=brown>$min_ratio\n";
print "<td>Max enrichment in repeats (ratio)<td><font color=brown>$max_repeat\n";
print "<tr><td>Match threshold for clustering<td><font color=brown>$match_thresh\n";
print "<td>Max number of motifs<td><font color=brown>$Nmatrix\n";
print "<tr><td>Clustering method<td><font color=brown>$method_list[$method]\n";
print "</table>\n";
if($loginname !~ /^guest/){
	print "(to change parameters use 'Edit project' button)\n";
}
print "<h3>Analysis Tools</h3>\n";
print "<TABLE border=0>\n";
print "<TR><TD><INPUT NAME=\"get_motif_button\" TYPE=\"button\"      VALUE=\"   Identify motifs      \" onClick=\"do_analysis('motif_get');\">\n";
print "<TD>Identify motifs over-repretented in sequence file #1 compared to sequence file #2 (#2 is optional)\n";
if($file_motif1){
	print "<TR><TD><INPUT NAME=\"cluster_motif_button\" TYPE=\"button\"  VALUE=\"   Cluster motifs      \" onClick=\"do_analysis('motif_cluster');\">\n";
	print "<TD>Cluster motifs from motif file based on similarity\n";
	print "<TR><TD><INPUT NAME=\"compare_motif_button\" TYPE=\"button\"  VALUE=\"   Compare motifs   \" onClick=\"do_analysis('motif_compare');\">\n";
	print "<TD>Compare motifs from motif file with other motif files (e.g. for annotation)\n";
	print "<TR><TD><INPUT NAME=\"explore_motif_button\" TYPE=\"button\"  VALUE=\"   Show motifs        \" onClick=\"explore_motifs();\">\n";
	print "<TD>Shows motif logo and PFM for motif file (includes motif copying)\n";
	print "<TR><TD><INPUT NAME=\"search_motif_button\" TYPE=\"button\"   VALUE=\"   Search motifs      \" onClick=\"do_analysis('motif_search');\">\n";
	print "<TD>Search sequence file #1 for occurence of motifs from motif file\n";
}
if($file_search){
	print "<TR><TD><INPUT NAME=\"search_results_button\" TYPE=\"button\" VALUE=\"Show search results\" onClick=\"do_analysis('show_search_results');\">\n";
	print "<TD>Show search results (motifs, sequences, frequency distribution)\n";
}
print "</TABLE>\n";
print "<HR NOSHADE></HR>\n";
#if($loginname !~ /^guest/){
#	print "<h3>Save project as a stand-alone web page</h3>\n";
#	print "File name: <INPUT NAME=\"web_page_name\" SIZE=20 VALUE=\"$dataname\"> (no special characters) &nbsp; &nbsp; &nbsp; &nbsp;\n";
#	print "<INPUT NAME=\"save_page_button\" TYPE=\"button\" VALUE=\"Save as web page\" onClick=\"save_page();\"><br>\n";
#}
print "</FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  project_edit
#**************************************
{
my $dataname = shift;
my @motif_list;
my @fasta_list;
my @repeat_list;
my @search_list;
my %hashEdit;
my %hashDefault;
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
my $found = 0;
while(my $line = <INFO>){
	my %hash=();
	if(length($line) < 3) { next; }
	read_config_line($line,\%hash);
	my $file_motif = $hash{"type_motif"};
	if($file_motif){ push(@motif_list,[$file_motif,$hash{"description"}]); next; }
	my $file_repeat = $hash{"type_repeat"};
	if($file_repeat){ push(@repeat_list,[$file_repeat,$hash{"description"}]); next; }
	my $file_fasta = $hash{"type_fasta"};
	if($file_fasta){ push(@fasta_list,[$file_fasta,$hash{"description"}]); next; }
	my $file_search = $hash{"type_search"};
	if($file_search){ push(@search_list,[$file_search,$hash{"description"}]); next; }
	my $data = $hash{"data"};
	if($hash{"mainPage"}){ %hashDefault = %hash; }
	elsif($data){
		if($dataname eq $data){
			$found=1;
			%hashEdit = %hash;
		}
	}
}
close INFO;
if(!$found){ %hashEdit=%hashDefault; }
@motif_list = sort {lc($a->[0]) cmp lc($b->[0])} @motif_list;
@repeat_list = sort {lc($a->[0]) cmp lc($b->[0])} @repeat_list;
@fasta_list = sort {lc($a->[0]) cmp lc($b->[0])} @fasta_list;
@search_list = sort {lc($a->[0]) cmp lc($b->[0])} @search_list;
if($loginname ne "public" && open(INFO, "$PATH_INFO/public-config.txt")){
	my @motif_list1=();
	my @fasta_list1=();
	my @repeat_list1=();
	my @search_list1=();
	while(my $line = <INFO>){
		my %hash=();
		read_config_line($line,\%hash);
		my $file_motif = $hash{"type_motif"};
		if($file_motif){ push(@motif_list1,["public-$file_motif",$hash{"description"}]); next; }
		my $file_repeat = $hash{"type_repeat"};
		if($file_repeat){ push(@repeat_list1,["public-$file_repeat",$hash{"description"}]); next; }
		my $file_fasta = $hash{"type_fasta"};
		if($file_fasta){ push(@fasta_list1,["public-$file_fasta",$hash{"description"}]); }
		my $file_search = $hash{"type_search"};
		if($file_search){ push(@search_list1,["public-$file_search",$hash{"description"}]); }
	}
	push(@motif_list, sort {lc($a->[0]) cmp lc($b->[0])} @motif_list1);
	push(@repeat_list, sort {lc($a->[0]) cmp lc($b->[0])} @repeat_list1);
	push(@fasta_list, sort {lc($a->[0]) cmp lc($b->[0])} @fasta_list1);
	push(@search_list, sort {lc($a->[0]) cmp lc($b->[0])} @search_list1);
}
close INFO;

my $file_test = $hashEdit{"file_test"};
my $file_control = $hashEdit{"file_control"};
my $file_repeat = $hashEdit{"file_repeat"};
my $file_motif1 = $hashEdit{"file_motif1"};
my $file_search = $hashEdit{"file_search"};

print "<HTML><HEAD><TITLE>CisFinder</TITLE>\n";
if($dataname && !$found){ error_message("Data $dataname not found"); }
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
my ($items,$descriptions) = get_array_lists(\@motif_list);
print "motif_list = new Array($items);\n";
print "motif_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@repeat_list);
print "repeat_list = new Array($items);\n";
print "repeat_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@fasta_list);
print "fasta_list = new Array($items);\n";
print "fasta_description = new Array($descriptions);\n";
($items,$descriptions) = get_array_lists(\@search_list);
print "search_list = new Array($items);\n";
print "search_description = new Array($descriptions);\n";
print "function update_description() {\n";
print "  var index = document.cisfinder.file_test.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_test.value = fasta_description[index-1]; }\n";
print "  else{ document.cisfinder.description_test.value = \"\"; }\n";
print "  index = document.cisfinder.file_control.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_control.value = fasta_description[index-1]; }\n";
print "  else{ document.cisfinder.description_control.value = \"\"; }\n";
print "  index = document.cisfinder.file_repeat.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_repeat.value = repeat_description[index-1]; }\n";
print "  else{ document.cisfinder.description_repeat.value = \"\"; }\n";
print "  index = document.cisfinder.file_motif1.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_motif1.value = motif_description[index-1]; }\n";
print "  else{ document.cisfinder.description_motif1.value = \"\"; }\n";
print "  var index = document.cisfinder.file_search.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_search.value = search_description[index-1]; }\n";
print "  else{ document.cisfinder.description_search.value = \"\"; }\n";
print "}\n";
print "function project_update() {\n";
print "	if(!document.cisfinder.file_test.value){ alert(\"You need to select a test file\"); return false; }\n";
print "	if(!document.cisfinder.data.value){ alert(\"You need to enter data name\"); return false; }\n";
print "	clear_marks();\n";
print "	document.cisfinder.action.value = \"show_project\";\n";
print "	document.cisfinder.update.value = \"show_project\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function cancel_edit() {\n";
print "	clear_marks();\n";
print "	document.cisfinder.action.value = \"continue\";\n";
print "	document.cisfinder.data.value = \"\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function edit_refresh() {\n";
print "	clear_marks();\n";
print "	document.cisfinder.action.value = \"project_edit\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function clear_marks() {\n";
print "	document.cisfinder.target = \"\"\n";
print "	document.cisfinder.file_download.value = \"\";\n";
print "	document.cisfinder.update.value = \"\";\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header("update_description();");
print "HELP: <a href=\"$HOME_ADDRESS/cisfinder-help.html\" TARGET=\"_blank324\">A guide to CisFinder</a><br>\n";

print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST LANGUAGE=\"javascript\" onsubmit=\"return cisfinder_onsubmit()\">\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"project_edit\">\n";
print "<INPUT NAME=\"file_download\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"update\" TYPE=\"hidden\" VALUE=\"\">\n";
my ($description,$use_repeats,$strand,$score,$FDR,$adjust_cg,$presence,$min_ratio,$max_repeat,$match_thresh,$method,$interval,$Nmatrix);
if($dataname){
	print "<p><H3>Edit project <font color=brown size=+1>$dataname</font></H3>\n";
	print "<INPUT NAME=\"update_conf\" TYPE=\"button\" VALUE=\"Update project\" onClick=\"project_update();\">\n";
	print "<INPUT NAME=\"update_abort\" TYPE=\"button\" VALUE=\"Cancel\" onClick=\"cancel_edit();\">\n";
	print "<INPUT NAME=\"button_edit_refresh\" TYPE=\"button\" VALUE=\"Refresh\" onClick=\"edit_refresh();\"><br>\n";
	print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
	$description = $hashEdit{"description"};
	$use_repeats = $hashEdit{"use_repeats"};
	$strand = $hashEdit{"strand"};
	$score = $hashEdit{"score"};
	$FDR = $hashEdit{"FDR"};
	$adjust_cg = $hashEdit{"adjust_cg"};
	$presence = $hashEdit{"presence"};
	$min_ratio = $hashEdit{"min_ratio"};
	$max_repeat = $hashEdit{"max_repeat"};
	$match_thresh = $hashEdit{"match_thresh"};
	$method = $hashEdit{"method"};
	$interval = $hashEdit{"interval"};
	$Nmatrix = $hashEdit{"Nmatrix"};
	if(!$Nmatrix){ $Nmatrix=500; }
}else{
	print "<p><H3>Create new project</H3>\n";
	print "<INPUT NAME=\"update_conf\" TYPE=\"button\" VALUE=\"Save project\" onClick=\"project_update();\">\n";
	print "<INPUT NAME=\"update_abort\" TYPE=\"button\" VALUE=\"Cancel\" onClick=\"cancel_edit();\">\n";
	print "<INPUT NAME=\"button_edit_refresh\" TYPE=\"button\" VALUE=\"Refresh\" onClick=\"edit_refresh();\"><br>\n";
	print "<INPUT NAME=\"data\" SIZE=40 VALUE=\"\"> Data name<br>\n";
	$score = 1;
	$FDR = 0.05;
	$match_thresh=0.7;
	$min_ratio=1.5;
	$max_repeat=1000;
	$Nmatrix = 500;
}
print "<INPUT NAME=\"description\" SIZE=40 VALUE=\"$description\"> Data description\n";
print "<p><b>Select files for analysis</b>\n";
print "<TABLE border=0>\n";
print "<TR><TD><b>File<TD><b>Description<TD><b>File type\n";
print "<TR><TD><select name=\"file_test\" onChange=update_description();>\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@fasta_list; ++$i){ 
	print "<option value=\"$fasta_list[$i]->[0]\"";
	if($fasta_list[$i]->[0] eq $file_test){ print " selected"; }
	print "> $fasta_list[$i]->[0]\n";
}
print "</select><td><INPUT NAME=\"description_test\" SIZE=40><TD>Sequence file #1 (test)\n";
print "<td><INPUT NAME=\"file_test_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_test');\">\n";
print "<TR><TD><select name=\"file_control\" onChange=update_description();>\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@fasta_list; ++$i){ 
	print "<option value=\"$fasta_list[$i]->[0]\"";
	if($fasta_list[$i]->[0] eq $file_control){ print " selected"; }
	print "> $fasta_list[$i]->[0]\n";
}
print "</select><td><INPUT NAME=\"description_control\" SIZE=40><TD> Sequence file #2 (control)\n";
print "<td><INPUT NAME=\"file_control_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_control');\">\n";
print "<TR><TD><select name=\"file_repeat\" onChange=update_description();>\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@repeat_list; ++$i){ 
	print "<option value=\"$repeat_list[$i]->[0]\"";
	if($repeat_list[$i]->[0] eq $file_repeat){ print " selected"; }
	print "> $repeat_list[$i]->[0]\n";
}
print "</select><td><INPUT NAME=\"description_repeat\" SIZE=40><TD>Repeat file\n";
print "<td><INPUT NAME=\"file_repeat_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_repeat');\">\n";
print "<TR><TD><select name=\"file_motif1\" onChange=update_description();>\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@motif_list; ++$i){ 
	print "<option value=\"$motif_list[$i]->[0]\"";
	if($motif_list[$i]->[0] eq $file_motif1){ print " selected"; }
	print "> $motif_list[$i]->[0]\n";
}
print "</select><td><INPUT NAME=\"description_motif1\" SIZE=40><TD>Motif file\n";
print "<td><INPUT NAME=\"file_motif_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_motif');\">\n";
print "<TR><TD><select name=\"file_search\" onChange=update_description();>\n";
print "<option value=\"\"> ------------------------------- select -------------------------------\n";
for(my $i=0; $i<@search_list; ++$i){ 
	print "<option value=\"$search_list[$i]->[0]\"";
	if($search_list[$i]->[0] eq $file_search){ print " selected"; }
	print "> $search_list[$i]->[0]\n";
}
print "</select><td><INPUT NAME=\"description_search\" SIZE=40><TD>Search results\n";
print "<td><INPUT NAME=\"file_search_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_search');\">\n";
print "</TABLE>\n";
print "<b>Parameters</b>\n";
my @yesno = ("No","Yes");
my @strand_list = ("Both","Forward","Forward vs. back");
my @score_list = ("Z","Z+ratio","Z+info","Z+ratio+info","Z+repeat","Z+ratio+repeat","Z+info+repeat","Z+ratio+info+repeat");
my @method_list = ("Similarity","Similarity and linkage","Linkage");
my @FDR_list = (1,0.5,0.2,0.1,0.05,0.01,0.001,0.0001);
my @min_ratio = (1,1.2,1.5,2,3,4,5);
my @max_repeat = (1,2,4,8,20,100,1000);
my @match_thresh = (0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.98);
my @interval_list = (5,10,20,50,100,200,500,1000,2000);
my @Nmatrix_list = (100,250,500,1000,2000,1000000);
print "<table border=0>\n";
print "<TR><TD><select name=\"use_repeats\">\n";
for(my $i=0; $i<2; ++$i){ 
	print "<option value=$i";
	if($use_repeats==$i){ print " selected"; }
	print "> $yesno[$i]\n";
}
print "</select><TD> Use repeats for search\n";
print "<TD><select name=\"strand\">\n";
for(my $i=0; $i<@strand_list; ++$i){ 
	print "<option value=$i";
	if($strand==$i){ print " selected"; }
	print "> $strand_list[$i]\n";
}
print "</select><TD> Search in strands\n";
print "<TR><TD><select name=\"FDR\">\n";
for(my $i=0; $i<@FDR_list; ++$i){ 
	print "<option value=$FDR_list[$i]";
	if($FDR==$FDR_list[$i]){ print " selected"; }
	print "> $FDR_list[$i]\n";
}
print "</select><TD> FDR\n";
print "<TD><select name=\"adjust_cg\">\n";
for(my $i=0; $i<2; ++$i){ 
	print "<option value=$i";
	if($adjust_cg==$i){ print " selected"; }
	print "> $yesno[$i]\n";
}
print "</select><TD> Adjust for CG/AT ratio and CpG\n";
print "<TR><TD><select name=\"presence\">\n";
for(my $i=0; $i<2; ++$i){ 
	print "<option value=$i";
	if($presence==$i){ print " selected"; }
	print "> $yesno[$i]\n";
}
print "</select><TD> Count motif once per sequence\n";
print "<TD><select name=\"score\">\n";
for(my $i=0; $i<@score_list; ++$i){ 
	print "<option value=$i";
	if($score==$i){ print " selected"; }
	print "> $score_list[$i]\n";
}
print "</select><TD> Score motifs by\n";
print "<TR><TD><select name=\"min_ratio\">\n";
for(my $i=0; $i<@min_ratio; ++$i){ 
	print "<option value=$min_ratio[$i]";
	if($min_ratio==$min_ratio[$i]){ print " selected"; }
	print "> $min_ratio[$i]\n";
}
print "</select><TD> Minimum enrichment ratio (test vs. control)\n";
print "<TD><select name=\"max_repeat\">\n";
for(my $i=0; $i<@max_repeat; ++$i){ 
	print "<option value=$max_repeat[$i]";
	if($max_repeat==$max_repeat[$i]){ print " selected"; }
	print "> $max_repeat[$i]\n";
}
print "</select><TD> Maximum enrichment in repeats (ratio)\n";
print "<TR><TD><select name=\"match_thresh\">\n";
for(my $i=0; $i<@match_thresh; ++$i){ 
	print "<option value=$match_thresh[$i]";
	if($match_thresh==$match_thresh[$i]){ print " selected"; }
	print "> $match_thresh[$i]\n";
}
print "</select><TD> Match threshold for clustering\n";
print "<TD><select name=Nmatrix>\n";
for(my $i=0; $i<@Nmatrix_list; ++$i){ 
	my $NNN = $Nmatrix_list[$i];
	print "<option value=$NNN";
	if($Nmatrix==$NNN){ print " selected"; }
	print "> $NNN\n";
}
print "</select><TD> Maximum number of motifs to find\n";
print "<TR><TD><select name=\"method\">\n";
for(my $i=0; $i<@method_list; ++$i){ 
	print "<option value=$i";
	if($method==$i){ print " selected"; }
	print "> $method_list[$i]\n";
}
print "</select><TD> Clustering method\n";
print "<TD><select name=\"interval\">\n";
for(my $i=0; $i<@interval_list; ++$i){ 
	print "<option value=$interval_list[$i]";
	if($interval==$interval_list[$i]){ print " selected"; }
	print "> $interval_list[$i]\n";
}
print "</select><TD> Interval (bp) for frequency distribution\n";
print "</table><p>\n";
print "</FORM>\n";
print "<HR NOSHADE></HR>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#************************************
sub  save_page
#************************************
{
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
my %hashProject;
my %fileDescription;
while(my $line = <INFO>){
	my %hash=();
	if(length($line) < 3) { next; }
	read_config_line($line,\%hash);
	if($hash{"mainPage"}){ next; }
	my $data = $hash{"data"};
	if($data eq $dataname){ %hashProject = %hash; last; }
}
close INFO;
my $file_test = $hashProject{"file_test"};
my $file_control = $hashProject{"file_control"};
my $file_repeat = $hashProject{"file_repeat"};
my $file_motif1 = $hashProject{"file_motif1"};
my $file_search = $hashProject{"file_search"};
my $description = $hashProject{"description"};
my $use_repeats = $hashProject{"use_repeats"};
my $FDR = $hashProject{"FDR"};
my $strand = $hashProject{"strand"};
my $score = $hashProject{"score"};
my $adjust_cg = $hashProject{"adjust_cg"};
my $presence = $hashProject{"presence"};
my $min_ratio = $hashProject{"min_ratio"};
my $max_repeat = $hashProject{"max_repeat"};
my $match_thresh = $hashProject{"match_thresh"};
my $method = $hashProject{"method"};
my $interval = $hashProject{"interval"};
my $Nmatrix = $hashProject{"Nmatrix"};
if(!$Nmatrix){ $Nmatrix=500; }

my ($description_test,$description_control,$description_repeat,$description_motif1,$description_search);
open(INFO, "$PATH_INFO/$loginname-config.txt");
while(my $line = <INFO>){
	my %hash=();
	if(length($line) < 3) { next; }
	read_config_line($line,\%hash);
	if($file_test && $file_test eq $hash{"type_fasta"}){ $description_test = $hash{"description"}; next; }
	if($file_control && $file_control eq $hash{"type_fasta"}){ $description_control = $hash{"description"}; next; }
	if($file_repeat && $file_repeat eq $hash{"type_repeat"}){ $description_repeat = $hash{"description"}; next; }
	if($file_motif1 && $file_motif1 eq $hash{"type_motif"}){ $description_motif1 = $hash{"description"}; next; }
	if($file_search && $file_search eq $hash{"type_search"}){ $description_search = $hash{"description"}; next; }
}
close INFO;
if($file_test !~ /^public-/){ $file_test = "$loginname-$file_test"; }
if($file_control && $file_control !~ /^public-/){ $file_control = "$loginname-$file_control"; }
if($file_repeat && $file_repeat !~ /^public-/){ $file_repeat = "$loginname-$file_repeat"; }
if($file_motif1 && $file_motif1 !~ /^public-/){ $file_motif1 = "$loginname-$file_motif1"; }
if($file_search && $file_search !~ /^public-/){ $file_search = "$loginname-$file_search"; }

my $web_page_name = $hashInput{"web_page_name"};
if(!$web_page_name){ error_message("No web page name"); }
if(!open(OUT, ">$PATH_HOME/saved/$loginname-$web_page_name.html")){ error_message("Cannt write to web page"); }
print OUT "<HTML><HEAD><TITLE>CisFinder</TITLE>\n";
print OUT "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print OUT "<!--\n";
print OUT "function clear_marks() {\n";
print OUT "	document.cisfinder.target = \"\"\n";
print OUT "	document.cisfinder.file_download.value = \"\";\n";
print OUT "}\n";
print OUT "function do_analysis(command) {\n";
print OUT "	clear_marks();\n";
print OUT "	var x = Math.round(Math.random()*1000);\n";
print OUT "	document.cisfinder.target = \"_BLANK\"+x;\n";
print OUT "	document.cisfinder.action.value = command;\n";
print OUT "	document.cisfinder.submit();\n";
print OUT "}\n";
print OUT "function download_file(file_type) {\n";
print OUT "	var filename;\n";
print OUT "	if(file_type == \"file_test\"){ filename = document.cisfinder.file_test.value; }\n";
print OUT "	if(file_type == \"file_control\"){ filename = document.cisfinder.file_control.value; }\n";
print OUT "	if(file_type == \"file_motif1\"){ filename = document.cisfinder.file_motif1.value; }\n";
print OUT "	if(file_type == \"file_repeat\"){ filename = document.cisfinder.file_repeat.value; }\n";
print OUT "	if(file_type == \"file_search\"){ filename = document.cisfinder.file_search.value; }\n";
print OUT "	if(!filename){ alert(\"File not selected\"); return false; }\n";
print OUT "	clear_marks();\n";
print OUT "	var x = Math.round(Math.random()*1000);\n";
print OUT "	document.cisfinder.target = \"_BLANK\"+x;\n";
print OUT "	document.cisfinder.action.value = \"\";\n";
print OUT "	document.cisfinder.file_download.value = filename;\n";
print OUT "	document.cisfinder.submit();\n";
print OUT "}\n";
print OUT "function data_manage(command){\n";
print OUT "	clear_marks();\n";
print OUT "	document.cisfinder.action.value = command;\n";
print OUT "	document.cisfinder.file_download.value = \"\";\n";
print OUT "	document.cisfinder.submit();\n";
print OUT "}\n";
print OUT "function explore_motifs() {\n";
print OUT "	clear_marks();\n";
print OUT "	var x = Math.round(Math.random()*1000);\n";
print OUT "	document.cisfinder.action.value = \"show_motifs\";\n";
print OUT "	document.cisfinder.target = \"_BLANK\"+x;\n";
print OUT "	document.cisfinder.submit();\n";
print OUT "}\n";
print OUT "// -->\n";
print OUT "</SCRIPT>\n";
my $text = get_header();
print OUT $text;
print OUT "HELP: <a href=\"$HOME_ADDRESS/cisfinder-help.html\" TARGET=\"_blank324\">A guide to CisFinder</a><br>\n";
print OUT "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST LANGUAGE=\"javascript\" onsubmit=\"return cisfinder_onsubmit()\">\n";
print OUT "<H3>Project: <font color=brown>$dataname</font></h3>\n";
print OUT "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"guest\">\n";
print OUT "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"analysis\">\n";
print OUT "<INPUT NAME=\"file_download\" TYPE=\"hidden\" VALUE=\"\">\n";
print OUT "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print OUT "<INPUT NAME=\"data_saved\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print OUT "<b>Description:</b> $description<p>\n";
print OUT "<TABLE border=1>\n";
print OUT "<tr><td><b>File type<td><b>File name<td><b>Description<td><b>Download\n";
print OUT "<tr><td><b>Sequence file #1 (test)<td><font color=brown>$file_test<td>$description_test\n";
print OUT "<INPUT NAME=\"file_test\" TYPE=\"hidden\" VALUE=\"$file_test\">\n";
print OUT "<td><INPUT NAME=\"file_test_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_test');\">\n";
if($file_control){
	print OUT "<tr><td><b>Sequence file #2 (control)<td><font color=brown>$file_control<td>$description_control\n";
	print OUT "<td><INPUT NAME=\"file_control_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_control');\">\n";
	print OUT "<INPUT NAME=\"file_control\" TYPE=\"hidden\" VALUE=\"$file_control\">\n";
}
if($file_repeat){
	print OUT "<tr><td><b>Repeat file<td><font color=brown>$file_repeat<td>$description_repeat\n";
	print OUT "<td><INPUT NAME=\"file_repeat_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_repeat');\">\n";
	print OUT "<INPUT NAME=\"file_repeat\" TYPE=\"hidden\" VALUE=\"$file_repeat\">\n";
}
if($file_motif1){
	print OUT "<tr><td><b>Motif file<td><font color=brown>$file_motif1<td>$description_motif1\n";
	print OUT "<td><INPUT NAME=\"file_motif1_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_motif1');\">\n";
	print OUT "<INPUT NAME=\"file_motif1\" TYPE=\"hidden\" VALUE=\"$file_motif1\">\n";
}
if($file_search){
	print OUT "<tr><td><b>Search results<td><font color=brown>$file_search<td>$description_search\n";
	print OUT "<td><INPUT NAME=\"file_search_download\" TYPE=\"button\" VALUE=\"Download\" onClick=\"download_file('file_search');\">\n";
	print OUT "<INPUT NAME=\"file_search\" TYPE=\"hidden\" VALUE=\"$file_search\">\n";
}
print OUT "</TABLE>\n";

print OUT "<h3>Parameters</h3>\n";
my @yesno = ("No","Yes");
my @strand_list = ("Both","Forward","Forward vs. back");
my @score_list = ("Z","Z+ratio","Z+info","Z+ratio+info","Z+repeat","Z+ratio+repeat","Z+info+repeat","Z+ratio+info+repeat");
my @method_list = ("Similarity","Similarity and linkage","Linkage");
print OUT "<table border=1>\n";
print OUT "<tr><td><b>Parameter<td><b>Value<td><b>Parameter<td><b>Value\n";
print OUT "<tr><td>Search repeats<td><font color=brown>$yesno[$use_repeats]\n";
print OUT "<INPUT NAME=\"use_repeats\" TYPE=\"hidden\" VALUE=\"$use_repeats\">\n";
print OUT "<td>Search in strands<td><font color=brown>$strand_list[$strand]\n";
print OUT "<INPUT NAME=\"strand\" TYPE=\"hidden\" VALUE=\"$strand\">\n";
print OUT "<tr><td>False Discovery Rate (FDR)<td><font color=brown>$FDR\n";
print OUT "<INPUT NAME=\"FDR\" TYPE=\"hidden\" VALUE=\"$FDR\">\n";
print OUT "<td>Adjust for CG/AT ratio and CpG<td><font color=brown>$yesno[$adjust_cg]\n";
print OUT "<INPUT NAME=\"adjust_cg\" TYPE=\"hidden\" VALUE=\"$adjust_cg\">\n";
print OUT "<tr><td>Count motif once per sequence<td><font color=brown>$yesno[$presence]\n";
print OUT "<INPUT NAME=\"presence\" TYPE=\"hidden\" VALUE=\"$presence\">\n";
print OUT "<td>Score motifs by<td><font color=brown>$score_list[$score]\n";
print OUT "<INPUT NAME=\"score\" TYPE=\"hidden\" VALUE=\"$score\">\n";
print OUT "<tr><td>Min enrichment ratio<td><font color=brown>$min_ratio\n";
print OUT "<INPUT NAME=\"min_ratio\" TYPE=\"hidden\" VALUE=\"$min_ratio\">\n";
print OUT "<td>Max enrichment in repeats (ratio)<td><font color=brown>$max_repeat\n";
print OUT "<INPUT NAME=\"max_repeat\" TYPE=\"hidden\" VALUE=\"$max_repeat\">\n";
print OUT "<tr><td>Match threshold for clustering<td><font color=brown>$match_thresh\n";
print OUT "<INPUT NAME=\"match_thresh\" TYPE=\"hidden\" VALUE=\"$match_thresh\">\n";
print OUT "<td>Max number of motifs<td><font color=brown>$Nmatrix\n";
print OUT "<INPUT NAME=\"Nmatrix\" TYPE=\"hidden\" VALUE=\"$Nmatrix\">\n";
print OUT "<tr><td>Clustering method<td><font color=brown>$method_list[$method]\n";
print OUT "<INPUT NAME=\"method\" TYPE=\"hidden\" VALUE=\"$method\">\n";
print OUT "</table>\n";
print OUT "<h3>Analysis Tools</h3>\n";
print OUT "<TABLE border=0>\n";
print OUT "<TR><TD><INPUT NAME=\"get_motif_button\" TYPE=\"button\"      VALUE=\"   Identify motifs      \" onClick=\"do_analysis('motif_get');\">\n";
print OUT "<TD>Identify motifs over-repretented in sequence file #1 compared to sequence file #2 (#2 is optional)\n";
if($file_motif1){
	print OUT "<TR><TD><INPUT NAME=\"cluster_motif_button\" TYPE=\"button\"  VALUE=\"   Cluster motifs      \" onClick=\"do_analysis('motif_cluster');\">\n";
	print OUT "<TD>Cluster motifs from motif file based on similarity\n";
	print OUT "<TR><TD><INPUT NAME=\"compare_motif_button\" TYPE=\"button\"  VALUE=\"   Compare motifs   \" onClick=\"do_analysis('motif_compare');\">\n";
	print OUT "<TD>Compare motifs from motif file with other motif files (e.g. for annotation)\n";
	print OUT "<TR><TD><INPUT NAME=\"explore_motif_button\" TYPE=\"button\"  VALUE=\"   Show motifs        \" onClick=\"explore_motifs();\">\n";
	print OUT "<TD>Shows motif logo and PFM for motif file (includes motif copying)\n";
	print OUT "<TR><TD><INPUT NAME=\"search_motif_button\" TYPE=\"button\"   VALUE=\"   Search motifs      \" onClick=\"do_analysis('motif_search');\">\n";
	print OUT "<TD>Search sequence file #1 for occurence of motifs from motif file\n";
}
if($file_search){
	print OUT "<TR><TD><INPUT NAME=\"search_results_button\" TYPE=\"button\" VALUE=\"Show search results\" onClick=\"do_analysis('show_search_results');\">\n";
	print OUT "<TD>Show search results (motifs, sequences, frequency distribution)\n";
}
print OUT "</TABLE>\n";
print OUT "<HR NOSHADE></HR>\n";
print OUT "</FORM>\n";
print OUT "</BODY>\n";
print OUT "</HTML>\n";
close(OUT);

#update configuration
open(OUT, ">$PATH_INFO/$loginname-config1.txt");
open(INFO, "$PATH_INFO/$loginname-config.txt");
while(my $line = <INFO>){
	if($line !~ /^saved_page=$web_page_name\t/){ print OUT $line; }
}
close INFO;
print OUT "saved_page=$web_page_name\tdescription=$description\n";
close OUT;
system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
system("rm $PATH_INFO/$loginname-config1.txt");

# Print page header
print "<HTML><HEAD><TITLE>CisFinder saved $web_page_name</TITLE>\n";
print_header();
print "<h3>Web page <a href=../saved/$loginname-$web_page_name.html>$web_page_name</a> is saved</h3>\n";
print "<INPUT NAME=\"close_button\" TYPE=button VALUE=\"Close window\" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  project_update
#**************************************
{
if(!$dataname){
	print "Warning: Dataname is empty in project_update";
	return;
}
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
if(!open(OUT, ">$PATH_INFO/$loginname-config1.txt")){
	error_message("Cannot update configuration file!");
}
while(my $line = <INFO>){
	my %hash=();
	if(length($line)<3){ next; }
	read_config_line($line,\%hash);
	if($hash{"mainPage"} || $dataname ne $hash{"data"}){ print OUT $line; }
}
close INFO;
my $file_test = $hashInput{"file_test"};
my $file_control = $hashInput{"file_control"};
my $file_repeat = $hashInput{"file_repeat"};
my $file_motif1 = $hashInput{"file_motif1"};
my $file_search = $hashInput{"file_search"};
my $description = $hashInput{"description"};
my $use_repeats = $hashInput{"use_repeats"};
my $strand = $hashInput{"strand"};
my $score = $hashInput{"score"};
my $FDR = $hashInput{"FDR"};
my $adjust_cg = $hashInput{"adjust_cg"};
my $presence = $hashInput{"presence"};
my $min_ratio = $hashInput{"min_ratio"};
my $max_repeat = $hashInput{"max_repeat"};
my $match_thresh = $hashInput{"match_thresh"};
my $method = $hashInput{"method"};
my $interval = $hashInput{"interval"};
my $Nmatrix = $hashInput{"Nmatrix"};
if(!$Nmatrix){ $Nmatrix=500; }
print OUT "data=$dataname\tfile_test=$file_test";
if($file_control){ print OUT "\tfile_control=$file_control"; }
if($file_repeat){ print OUT "\tfile_repeat=$file_repeat"; }
if($file_motif1){ print OUT "\tfile_motif1=$file_motif1"; }
if($file_search){ print OUT "\tfile_search=$file_search"; }
print OUT "\tdescription=$description\tuse_repeats=$use_repeats\tstrand=$strand";
print OUT "\tscore=$score\tFDR=$FDR\tadjust_cg=$adjust_cg\tpresence=$presence\tmin_ratio=$min_ratio";
print OUT "\tmax_repeat=$max_repeat\tmatch_thresh=$match_thresh\tmethod=$method\tinterval=$interval\tNmatrix=$Nmatrix\n";
close OUT;
system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
system("rm $PATH_INFO/$loginname-config1.txt");
return;
}

#**************************************
sub  check_configuration
#**************************************
{
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
if(!open(OUT, ">$PATH_INFO/$loginname-config1.txt")){
	error_message("Cannot update configuration file!");
}
my %hashDefault=();
my $changed=0;
while(my $line = <INFO>){
	if(length($line)<3){ next; }
	my %hash=();
	read_config_line($line,\%hash);
	if($hash{"mainPage"}){
		%hashDefault = %hash;
	}else{
		print OUT $line;
	}
}
close INFO;
foreach my $key (keys %hashInput){
	if($key eq "data" || $key =~/mainPage|action|loginname|passwd|saved_page|description|delete|download|submit|select_|fasta_file|upload_|motif_num|load_file|pasted_/){ next; }
	my $value = $hashDefault{$key};
	if($value ne $hashInput{$key}){
		$hashDefault{$key} = $hashInput{$key};
		$changed=1;
	}
}
if($changed){
	print OUT "mainPage=mainPage";
	foreach my $key (keys %hashDefault){
		if($key eq "data" || $key =~/mainPage|action|loginname|passwd|saved_page|description|delete|download|submit|select_|fasta_file|upload_|motif_num|load_file|pasted_/){ next; }
		print OUT "\t$key=$hashDefault{$key}";
	}
	print OUT "\n";
}
close OUT;
if($changed){
	system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
}
system("rm $PATH_INFO/$loginname-config1.txt");
return;
}

#**************************************
sub  project_delete
#**************************************
{
if(!$dataname){ error_message("Dataname is empty in project_delete"); }
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
if(!open(OUT, ">$PATH_INFO/$loginname-config1.txt")){
	error_message("Cannot update configuration file!");
}
while(my $line = <INFO>){
	my %hash=();
	if(length($line)<3){ next; }
	read_config_line($line,\%hash);
	if($hash{"mainPage"} || $dataname ne $hash{"data"}){ print OUT $line; }
}
close INFO;
close OUT;
system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
system("rm $PATH_INFO/$loginname-config1.txt");
return;
}

#**************************************
sub  file_delete
#**************************************
{
my $filename = $hashInput{"file_delete"};
my $saved_page = 0;
if($filename =~ /^saved\//){
	system("rm ../$filename.html");
	$saved_page = 1;
	$filename =~ s/^saved\/$loginname-//;
}else{
	system("rm $PATH_DATA/$loginname-$filename");
}
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
if(!open(OUT, ">$PATH_INFO/$loginname-config1.txt")){
	error_message("Cannot update configuration file!");
}
while(my $line = <INFO>){
	my %hash=();
	if(length($line)<3){ next; }
	read_config_line($line,\%hash);
	if($saved_page && $hash{"saved_page"} eq $filename){ next; }
	if($filename eq $hash{"type_fasta"}){
		my $file = $hash{"coordinates"};
		if($file){ system("rm $PATH_DATA/$loginname-$file"); }
		$file = $hash{"attributes"};
		if($file){ system("rm $PATH_DATA/$loginname-$file"); }
		$file = $hash{"conservation"};
		if($file){ system("rm $PATH_DATA/$loginname-$file"); }
	}
	if($filename eq $hash{"type_fasta"} || $filename eq $hash{"type_repeat"} || $filename eq $hash{"type_motif"} || $filename eq $hash{"type_search"}){ next; }
	if(!$hash{"mainPage"} && $hash{"data"} && $line =~ /=$filename(\t|$)/){
		$line =~ s/\t[^\t]+=$filename\t/\t/g;
		$line =~ s/\t[^\t]+=$filename\n/\n/;
	}
	print OUT $line;
}
close INFO;
close OUT;
system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
system("rm $PATH_INFO/$loginname-config1.txt");
return;
}

#**************************************
sub  file_public
#**************************************
{
my $filename = $hashInput{"file_public"};
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
my %hashNew;
my $file_type;
while(my $line = <INFO>){
	my %hash=();
	if(length($line)<3){ next; }
	read_config_line($line,\%hash);
	if($filename eq $hash{"type_fasta"}){
		%hashNew = %hash; $file_type = "type_fasta"; last;
	}elsif($filename eq $hash{"type_motif"}){
		%hashNew = %hash; $file_type = "type_motif"; last;
	}elsif($filename eq $hash{"type_repeat"}){
		%hashNew = %hash; $file_type = "type_repeat"; last;
	}elsif($filename eq $hash{"type_search"}){
		%hashNew = %hash; $file_type = "type_search"; last;
	}
}
close INFO;
if(!open(OUT, ">>$PATH_INFO/public-config.txt")){
	error_message("Cannot update configuration file!");
}
system("cp $PATH_DATA/$loginname-$filename $PATH_DATA/public-$filename");
print OUT "$file_type=$filename";
foreach my $key (keys %hashNew){
	if($key eq $file_type){ next; }
	print OUT "\t$key=$hashNew{$key}";
	if($file_type eq "type_fasta"){
		my $file = $hashNew{"coordinates"};
		if($file){ system("cp $PATH_DATA/$loginname-$file $PATH_DATA/public-$file"); }
		$file = $hashNew{"attributes"};
		if($file){ system("cp $PATH_DATA/$loginname-$file $PATH_DATA/public-$file"); }
		$file = $hashNew{"conservation"};
		if($file){ system("cp $PATH_DATA/$loginname-$file $PATH_DATA/public-$file"); }
	}
}
print OUT "\n";
close OUT;
return;
}

#**************************************
sub  file_download
#**************************************
{
my $filename = $hashInput{"file_download"};
if(!$hashInput{"data_saved"} && $filename !~ /^public-/){
	$filename = "$loginname-$filename";
}
my $outputID = get_outputID(1);
if(!open(INFO, "$PATH_DATA/$filename")){
	error_message("File $filename not found"); 
}else{
	close INFO;
}
`cp $PATH_DATA/$filename $PATH_OUTPUT/$outputID.txt`;
print_header();
print "<h2>Fil $filename is avalable below</h2>";
print "<a href=\"$HOME_ADDRESS/output/$outputID.txt\">$outputID.txt</a><p>";
print "Hint: use right-button mouse click to download<p>";
exit(0);
}

#**************************************
sub  motif_get
#**************************************
{
my %hashConfig=%hashInput;
if(!$hashInput{"data_saved"}){
	if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){ error_message("Profile not found"); }
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash=();
		read_config_line($line,\%hash);
		if($hashInput{"mainPage"}){
			if($hash{"mainPage"}){ %hashConfig = %hash; last; }
		}else{
			if($dataname && $dataname eq $hash{"data"}){ %hashConfig = %hash; last; }
		}
	}
	close INFO;
	if($hashInput{"mainPage"}){ $dataname=""; }
}
my $file_test = $hashConfig{"file_test"};
my $file_control = $hashConfig{"file_control"};
my $file_repeat = $hashConfig{"file_repeat"};
my $use_repeats = $hashConfig{"use_repeats"};
my $strand = $hashConfig{"strand"};
my $presence = $hashConfig{"presence"};
my $FDR = $hashConfig{"FDR"};
my $min_ratio = $hashConfig{"min_ratio"};
my $max_repeat = $hashConfig{"max_repeat"};
my $method = $hashConfig{"method"};
my $adjust_cg = $hashConfig{"adjust_cg"};
my $score = $hashConfig{"score"};
my $match_thresh = $hashConfig{"match_thresh"};
my $Nmatrix = $hashConfig{"Nmatrix"};
if(!$FDR){ $FDR=0.05; }
if(!$min_ratio){ $min_ratio=1.5; }
if(!$max_repeat){ $max_repeat=1000; }
if(!defined($score)){ $score=1; }
if(!$match_thresh){ $match_thresh=0.75; }
if(!$Nmatrix){ $Nmatrix=500; }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();

my $taskID = get_outputID(4);
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"motif_get1\">\n";
print "<INPUT NAME=\"taskID\" TYPE=\"hidden\" VALUE=\"$taskID\">\n";
#if($hashInput{"data_saved"}){
#	print "<INPUT NAME=\"data_saved\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
#}
if($hashInput{"mainPage"}){
	print "<INPUT NAME=\"mainPage\" TYPE=\"hidden\" VALUE=\"mainPage\">\n";
}
print "<INPUT NAME=\"file_test\" TYPE=\"hidden\" VALUE=\"$file_test\">\n";
print "<h3>Files selected for motif identification</h3>\n";
print "<TABLE border=1>\n";
print "<tr><td><b>File type<td><b>File name\n";
print "<tr><td><b>Sequence file #1 (test)<td><font color=brown>$file_test\n";
if($file_control){
	print "<tr><td><b>Sequence file #2 (control)<td><font color=brown>$file_control\n";
	print "<INPUT NAME=\"file_control\" TYPE=\"hidden\" VALUE=\"$file_control\">\n";
}
if($file_repeat){
	print "<tr><td><b>Repeat file<td><font color=brown>$file_repeat\n";
	print "<INPUT NAME=\"file_repeat\" TYPE=\"hidden\" VALUE=\"$file_repeat\">\n";
}
print "</TABLE>\n";

print "<h3>Parameters for motif identification (modify if needed)</h3>\n";
my @yesno = ("No","Yes");
my @strand_list = ("Both","Forward","Forward vs. back");
my @score_list = ("Z","Z+ratio","Z+info","Z+ratio+info","Z+repeat","Z+ratio+repeat","Z+info+repeat","Z+ratio+info+repeat");
my @method_list = ("Similarity","Similarity and linkage","Linkage");
my @FDR_list = (1,0.5,0.2,0.1,0.05,0.01,0.001,0.0001);
my @min_ratio = (1.2,1.5,2,3,4,5);
my @max_repeat = (1,2,4,8,20,100,1000);
my @match_thresh = (0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.98);
my @Nmatrix_list = (100,250,500,1000,2000,1000000);
print "<table border=0>\n";
print "<TR><TD><select name=\"use_repeats\">\n";
for(my $i=0; $i<2; ++$i){ 
	print "<option value=$i";
	if($use_repeats==$i){ print " selected"; }
	print "> $yesno[$i]\n";
}
print "</select><TD> Use repeats for search\n";
print "<TD><select name=\"strand\">\n";
for(my $i=0; $i<@strand_list; ++$i){ 
	print "<option value=$i";
	if($strand==$i){ print " selected"; }
	print "> $strand_list[$i]\n";
}
print "</select><TD> Search in strands\n";
print "<TR><TD><select name=\"FDR\">\n";
for(my $i=0; $i<@FDR_list; ++$i){ 
	print "<option value=$FDR_list[$i]";
	if($FDR==$FDR_list[$i]){ print " selected"; }
	print "> $FDR_list[$i]\n";
}
print "</select><TD> FDR\n";
print "<TD><select name=\"adjust_cg\">\n";
for(my $i=0; $i<2; ++$i){ 
	print "<option value=$i";
	if($adjust_cg==$i){ print " selected"; }
	print "> $yesno[$i]\n";
}
print "</select><TD> Adjust for CG/AT ratio and CpG\n";
print "<TR><TD><select name=\"presence\">\n";
for(my $i=0; $i<2; ++$i){ 
	print "<option value=$i";
	if($presence==$i){ print " selected"; }
	print "> $yesno[$i]\n";
}
print "</select><TD> Count motif once per sequence\n";
print "<TD><select name=\"score\">\n";
for(my $i=0; $i<@score_list; ++$i){ 
	print "<option value=$i";
	if($score==$i){ print " selected"; }
	print "> $score_list[$i]\n";
}
print "</select><TD> Score motifs by\n";
print "<TR><TD><select name=\"min_ratio\">\n";
for(my $i=0; $i<@min_ratio; ++$i){ 
	print "<option value=$min_ratio[$i]";
	if($min_ratio==$min_ratio[$i]){ print " selected"; }
	print "> $min_ratio[$i]\n";
}
print "</select><TD> Minimum enrichment ratio (test vs. control)\n";
print "<TD><select name=\"max_repeat\">\n";
for(my $i=0; $i<@max_repeat; ++$i){ 
	print "<option value=$max_repeat[$i]";
	if($max_repeat==$max_repeat[$i]){ print " selected"; }
	print "> $max_repeat[$i]\n";
}
print "</select><TD> Maximum enrichment in repeats (ratio)\n";
print "<TR><TD><select name=\"match_thresh\">\n";
for(my $i=0; $i<@match_thresh; ++$i){ 
	print "<option value=$match_thresh[$i]";
	if($match_thresh==$match_thresh[$i]){ print " selected"; }
	print "> $match_thresh[$i]\n";
}
print "</select><TD> Match threshold for clustering\n";
print "<TD><select name=Nmatrix>\n";
for(my $i=0; $i<@Nmatrix_list; ++$i){ 
	my $NNN = $Nmatrix_list[$i];
	print "<option value=$NNN";
	if($Nmatrix==$NNN){ print " selected"; }
	print "> $NNN\n";
}
print "</select><TD> Maximum number of motifs to find\n";
print "<TR><TD><select name=\"method\">\n";
for(my $i=0; $i<@method_list; ++$i){ 
	print "<option value=$i";
	if($method==$i){ print " selected"; }
	print "> $method_list[$i]\n";
}
print "</select><TD> Clustering method\n";
print "</table><p>\n";
print "<TR><TD><INPUT NAME=\"submit\" TYPE=submit VALUE=Continue>\n";
print "</FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  motif_get1
#**************************************
{
my $file_test = $hashInput{"file_test"};
my $file_control = $hashInput{"file_control"};
my $file_repeat = $hashInput{"file_repeat"};
my $use_repeats = $hashInput{"use_repeats"};
my $strand = $hashInput{"strand"};
my $score = $hashInput{"score"};
my $FDR = $hashInput{"FDR"};
my $adjust_cg = $hashInput{"adjust_cg"};
my $presence = $hashInput{"presence"};
my $min_ratio = $hashInput{"min_ratio"};
my $max_repeat = $hashInput{"max_repeat"};
my $method = $hashInput{"method"};
my $match_thresh = $hashInput{"match_thresh"};
my $taskID = $hashInput{"taskID"};
my $Nmatrix = $hashInput{"Nmatrix"};
if(!$Nmatrix){ $Nmatrix=500; }
if($file_test !~ /^public-/ && !$hashInput{"data_saved"}){ $file_test = "$loginname-$file_test"; }
if($file_control && $file_control !~ /^public-/ && !$hashInput{"data_saved"}){ $file_control = "$loginname-$file_control"; }
if($file_repeat && $file_repeat !~ /^public-/ && !$hashInput{"data_saved"}){ $file_repeat = "$loginname-$file_repeat"; }

my @motif_list=();
if(open(INFO, "$PATH_INFO/$loginname-config.txt")){
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash1=();
		read_config_line($line,\%hash1);
		if($hash1{"type_motif"}){ push(@motif_list,$hash1{"type_motif"}); }
	}
	close INFO;
	@motif_list = sort @motif_list;
}
my $list1 = join("\",\"",@motif_list);
if($list1){ $list1 = "\"".$list1."\""; }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "motif_list = new Array($list1);\n";
print "function save_motif_file(){\n";
print "	if(!document.cisfinder.file_motif.value){\n";
print "		alert(\"Enter name for new motif file\"); return(false);\n";
print "	}\n";
print "	var filename = document.cisfinder.file_motif.value;\n";
print "	if(filename.search(/ |\\/|\\&|\\t|\\%|\\+|\\@|\\!/) >= 0){\n";
print "		alert(\"File name should have neither spaces nor special characters\");\n";
print "		return(false);\n";
print "	}\n";
print "	var i;\n";
print "	for(i=0; i<motif_list.length; ++i){\n";
print "		if(filename == motif_list[i] && !confirm(\"Do you want to overwrite existing file \"+filename+\"?\")){\n";
print "			alert(\"Cancelled\");\n";
print "			return(false);\n";
print "		}\n";
print "	}\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.file_motif_output.value=document.cisfinder.file_motif_select.options[document.cisfinder.file_motif_select.selectedIndex].value;\n";
print "	document.cisfinder.action.value=\"save_motifs\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function show_motifs(file,motif_elem){\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value=\"show_motifs\";\n";
print "	document.cisfinder.file_motif_output.value=file;\n";
print "	document.cisfinder.file_motif_elem.value=motif_elem;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();

print "<b>Analyzing sequences...</b><br>\n";
print "If connection is lost, use Refresh button<br>\n";
my $command1="./get_motifs.pl $loginname $sessionID $file_test $taskID $FDR $min_ratio $Nmatrix $match_thresh $max_repeat";
if($file_control){ $command1 .= " -c $file_control"; }
if($file_repeat){  $command1 .= " -rep $file_repeat"; }
if($use_repeats){  $command1 .= " -userep"; }
if($strand){       $command1 .= " -strand $strand"; }
if($score){        $command1 .= " -score $score"; }
if($adjust_cg){    $command1 .= " -cg"; }
if($method){       $command1 .= " -pos $method"; }
if($presence){     $command1 .= " -one"; }

my $response="";
if(open (INFO, "$PATH_OUTPUT/$taskID.txt")){
	while(my $line=<INFO>){ $response .= $line; }
	close INFO;
}else{
	#print "$command1<p>\n";
	`$command1`;
	$response="";
	if(open (INFO, "$PATH_OUTPUT/$taskID.txt")){
		while(my $line=<INFO>){ $response .= $line; }
		close INFO;
	}
}
#print "RESPONSE: $response<p>\n";
if($response =~ /\<\/HTML\>/i){
	print $response;
}else{
	print "<b>Not finished yet... Try to refresh the page later.</b>\n";
}
exit(0);
}

#**************************************
sub  show_cluster
#**************************************
{
my $cluster = $hashInput{"motif"};
my $seqLogo = $hashInput{"seqLogo"};
my $outputFind = $hashInput{"file_motif_elem"};
my $file_motif_output = $hashInput{"file_motif_output"};

# Print page header
print "<HTML><HEAD><TITLE>Motif $cluster</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
my @motif_list=();
if(open(INFO, "$PATH_INFO/$loginname-config.txt")){
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash=();
		read_config_line($line,\%hash);
		if($hash{"type_motif"}){ push(@motif_list,$hash{"type_motif"}); }
	}
	close INFO;
}
my $list1 = join("\",\"",@motif_list);
if($list1){ $list1 = "\"".$list1."\""; }
print "motif_list = new Array($list1);\n";
print "function motif_check(){\n";
print "	if(!document.cisfinder.pfm_name.value){\n";
print "		alert(\"Enter name for new motif\"); return(false);\n";
print "	}\n";
print "	if(document.cisfinder.pfm_end.value-document.cisfinder.pfm_start.value < 6){\n";
print "		alert(\"Selected motif length is too short (<7)\"); return(false);\n";
print "	}\n";
print "	if(document.cisfinder.file_motif_save.selectedIndex == 0){\n";
print "		var filename = prompt(\"Enter file name for new motif file\",\"\");\n";
print "		if(!filename){ alert(\"Operation cancelled\"); return(false); }\n";
print "		var i;\n";
print "		for(i=0; i<motif_list.length; ++i){\n";
print "			if(filename == motif_list[i]){\n";
print "				alert(\"File name already exist!\"); return false;\n";
print "			}\n";
print "		}\n";
print "		var y=document.createElement('option');\n";
print "		y.text=filename;\n";
print "		y.value=filename;\n";
print "  	try{ document.cisfinder.file_motif_save.add(y,null); }\n";
print "  	catch(ex){ document.cisfinder.file_motif_save.add(y); }\n";
print "  	document.cisfinder.file_motif_save.selectedIndex = document.cisfinder.file_motif_save.length-1;\n";
print "	}\n";
print "	return true;\n";
print "}\n";
print "function motif_save(){\n";
print "	document.cisfinder.action.value = \"save_PFM\";\n";
print "	if(!motif_check()){ return false; }\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function motif_generate(){\n";
print "	document.cisfinder.action.value = \"generate_PFM\";\n";
print "	if(!motif_check()){ return false; }\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<h3>Motif $cluster</h3>\n";

if(!open(INFO, "$PATH_OUTPUT/$file_motif_output")){
	print "Cluster output file not found!\n"; exit(0);
}
my $line="";
my @headers;
my @headerInd;
my $junk;
while($line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^Headers/){
		($junk,$junk,@headers) = split("\t",$line);
	}	
	if($line =~ /^>$cluster\t/){ last; }
}
for(my $j=0; $j<@headersGlob; ++$j){
	$headerInd[$j]=-1;
	for(my $i=0; $i<@headers; ++$i){
		if($headersGlob[$j] eq $headers[$i]){ $headerInd[$j]=$i; }
	}
}
$line =~ s/^\>//;
my ($clName,@items) = split(/\t/,$line);
if($headerInd[2]>=0 && $items[$headerInd[2]]==1){
	my $i=0;
	while($i<@headers && $headers[$i] ne "MotifName"){ ++$i; }
	my $motifName = $items[$i];
	print "<b>Motif name: $motifName</b><p>\n";
}
my $pattern = $items[$headerInd[0]];
print "<TABLE BORDER=1><TR>\n";
for(my $j=1; $j<@headersGlob; ++$j){
	if($headerInd[$j]>=0){ print "<TD width=80><center>$headersGlob[$j]"; }
}
print "<TR>\n";
for(my $j=1; $j<@headersGlob; ++$j){
	if($headerInd[$j]<0){ next; }
	print "<TD><center>$items[$headerInd[$j]]\n";
}
print "</TABLE>\n";

print "<img src=../output/$seqLogo.gif border=0></a>\n";
print "<TABLE BORDER=0><TR>\n";
my(@ch) = split(//,$pattern);
my $len = @ch;
for(my $i=0; $i<@ch; ++$i){ print "<TD width=23><center>$ch[$i]\n"; }
print "<TR>\n";
for(my $i=1; $i<=@ch; ++$i){ print "<TD><center>$i\n"; }
print "</TABLE>\n";

print "<h3>Position Frequency Matrix (PFM)</h3>\n";
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<INPUT NAME=\"cluster\" TYPE=\"hidden\" VALUE=\"$cluster\">\n";
print "<INPUT NAME=\"motif\" TYPE=\"hidden\" VALUE=\"$cluster\">\n";
print "<INPUT NAME=\"outputFind\" TYPE=\"hidden\" VALUE=$outputFind>\n";
print "<INPUT NAME=\"file_motif_output\" TYPE=\"hidden\" VALUE=$file_motif_output>\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"save_PFM\">\n";

print "<TABLE BORDER=0>\n";
print "<TR><TD><TD><center>Start position<TD><center>End position<TD><center>Reverse\n";
print "<TD><center>Name<TD><center>Add PFM to file\n";
print "<TR><TD>\n";
print "<INPUT NAME=\"save_pfm_button\" TYPE=\"button\" VALUE=\"Save motif\" onClick=\"motif_save();\">";
print "<TD><center><INPUT NAME=\"pfm_start\" SIZE=2 VALUE=1>\n";
print "<TD><center><INPUT NAME=\"pfm_end\" SIZE=2 VALUE=$len>\n";
print "<TD><center><INPUT NAME=\"reverse_PFM\" TYPE=\"checkbox\">\n";
print "<TD><center><INPUT NAME=\"pfm_name\" SIZE=5>\n";
print "<TD><center><SELECT name=\"file_motif_save\">\n";
print "<option value=\"\"> New file\n";
for(my $i=0; $i<@motif_list; ++$i){ 
	print "<option value=\"$motif_list[$i]\"";
	if($i == @motif_list-1){ print " selected"; }
	print "> $motif_list[$i]\n";
}
print "</select>\n";
print "</TABLE>\n";

print "<TABLE BORDER=1><TR><TD><b><font color=blue>Position<TD width=70><center><b><font color=blue>A<TD width=70><center><b><font color=blue>C<TD width=70><center><b><font color=blue>G<TD width=70><center><b><font color=blue>T\n";
while($line = <INFO>){
	$line =~ s/\n$//;
	if(!$line){ last; }
	my @items = split("\t",$line);
	$items[0]++;
	print "<TR>";
	for(my $i=0; $i<5; ++$i){
		print "<TD><center>$items[$i]\n";
	}
}
print "</TABLE>\n";
while($line = <INFO>){
	if($line =~ /^>/ || $line =~ /^members/){ last; }
}
if(!$line || $line =~ /^>/){ exit(0); }

print "<HR NOSHADE COLOR=BLUE></HR>\n";
print "<h3>Member Motifs</h3>\n";
print "<INPUT NAME=\"generate_pfm_button\" TYPE=\"button\" VALUE=\"Generate PFM\" onClick=\"motif_generate();\">\n";
print "Use this button to generate PFM from motifs selected below<br>\n";

print "<TABLE BORDER=0><TR><TD><TD><b><font color=blue>Motif Name<TD><center><b><font color=blue>Pattern\n";
if($headerInd[3]>=0){ print "<TD width=80><center><b><font color=blue>Freq"; }
if($headerInd[4]>=0){ print "<TD width=80><center><b><font color=blue>Ratio"; }
if($headerInd[8]>=0){ print "<TD><b><font color=blue>Repeat ratio<TD><center><b><font color=blue>Repeat name\n"; }
print "\n";
while($line = <INFO>){
	$line =~ s/\n$//;
	if(!$line){ last; }
	my ($name,$pattern,$freq,$ratio,$dir,$offset,$repeat,$repeatName) = split("\t",$line);
	print "<TR><TD>";
	print "<INPUT NAME=\"use_motif_$name\" TYPE=\"checkbox\" CHECKED>";
	print "<TD>$name<TD><pre>$pattern";
	if($headerInd[3]>=0){ print "<TD><center>$freq"; }
	if($headerInd[4]>=0){ print "<TD><center>$ratio"; }
	if($headerInd[8]>=0 && $repeat){ print "<TD><center>$repeat"; }
	if($repeatName){ print "<TD><center>$repeatName"; }
	print "\n";
}
print "</TABLE>\n";
close INFO;
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  save_PFM
#**************************************
{
my $motif = $hashInput{"motif"};
my $pfm_start = $hashInput{"pfm_start"};
my $pfm_end = $hashInput{"pfm_end"};
my $file_motif_save = $hashInput{"file_motif_save"};
my $reverse_PFM = $hashInput{"reverse_PFM"};
my $pfm_name = $hashInput{"pfm_name"};
my $file_motif = $hashInput{"file_motif"};
my $file_motif_output = $hashInput{"file_motif_output"};

my $filenameMotif = "$PATH_DATA/$loginname-$file_motif";
my $file_motif_source = $file_motif;
if($file_motif =~ /^public-/){ $filenameMotif = "$PATH_DATA/$file_motif"; }
if($file_motif_output){
	$filenameMotif = "$PATH_OUTPUT/$file_motif_output";
	$file_motif_source = $file_motif_output;
}

# Print page header
print "<HTML><HEAD><TITLE>Motif $motif</TITLE>\n";
print_header();

if(!$file_motif_save){
	print "File motif to save is empty!\n"; exit(0);
}
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
if(!open(OUT, ">$PATH_INFO/$loginname-config1.txt")){
	print "Cannot update configuration file!\n"; exit(0);
}
my $found=0;
while(my $line = <INFO>){
	my %hash=();
	if(length($line)<3){ next; }
	read_config_line($line,\%hash);
	my $file_motif1 = $hash{"type_motif"};
	if($file_motif1 eq $file_motif_save){
		$found=1;
	}
	print OUT $line;
}
if(!$found){
	print OUT "type_motif=$file_motif_save\tdescription=Motifs from $file_motif_source\n";
}
close INFO;
close OUT;
if(!$found){
	system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
}
system("rm $PATH_INFO/$loginname-config1.txt");

if(!open(INFO, $filenameMotif)){
	print "Motif output file not found!\n"; exit(0);
}
if($found){
	open(INFO1, "$PATH_DATA/$loginname-$file_motif_save");
	my $last_line;
	while(my $line = <INFO1>){
		if($line =~ /^>$pfm_name\t/){
			print "Motif \'$pfm_name\' is present in file \'$file_motif_save\'! Use different name\n"; exit(0);
		}
		$last_line = $line;
	}
	close INFO1;
	open(OUT, ">>$PATH_DATA/$loginname-$file_motif_save");
	if(length($last_line) > 2){
		print OUT "\n";
	}
}else{
	open(OUT, ">$PATH_DATA/$loginname-$file_motif_save");
}
my @headers;
my $junk;
my $line;
my @headerInd;
while($line = <INFO>){
	if(!$found && $line =~ /^Headers|^Parameters/){
		print OUT $line;
		if($line =~ /^Headers/){
			$line =~ s/\n$//;
			($junk,@headers) = split("\t",$line);
		}
	}
	if($line =~ /^>$motif\t/){ last; }
}
if($line !~ /^>$motif\t/){ print "Motif $motif not found"; close OUT; exit(0); }
$line =~ s/^>//;
$line =~ s/\n$//;
my ($name,$pattern,$patternRev,$freq,$ratio,$info,$score,$FDR,$palindrome,$Nmembers,$method)=split(/\t/,$line);
for(my $j=0; $j<@headersGlob; ++$j){
	$headerInd[$j]=-1;
	for(my $i=0; $i<@headers; ++$i){
		if($headersGlob[$j] eq $headers[$i]){ $headerInd[$j]=$i; }
		#print "$headersGlob[$j] $headers[$i] $headerInd[$j]<br>\n";
	}
}
my @Matrix=();
my $count=0;
while($line = <INFO>){
	$line =~ s/\n$//;
	if(!$line){ last; }
	my ($pos,@m) = split(/\t/,$line);
	if($pos != $count){ print "Error in matrix format!\n"; exit(0); }
	push(@Matrix,\@m);
	++$count;
}
close INFO;
my $len = @Matrix;
my $lenNew = $pfm_end-$pfm_start+1;
$pattern = substr($pattern,$pfm_start-1,$lenNew);
$patternRev = substr($patternRev,$len-$pfm_end,$lenNew);

if($pfm_start > 1){
	splice(@Matrix,0,$pfm_start-1);
}
if(@Matrix > $lenNew){
	splice(@Matrix,$lenNew);
}
if($reverse_PFM){
	@Matrix = reverse(@Matrix);
	for(my $i=0; $i<$len; ++$i){ @{$Matrix[$i]} = reverse(@{$Matrix[$i]}); }
	($pattern,$patternRev) = ($patternRev,$pattern);
}
$info=0;
my $maxEntropy = 2;
my $log2 = log(2);
for(my $i=0; $i<$lenNew; ++$i){
	my $sum=0;
	for(my $k=0; $k<4; ++$k){
		$sum += $Matrix[$i]->[$k];
	}
	my $x=0;
	my $coef=1;
	my $entropy1=0;
	for(my $k=0; $k<4; ++$k){
		my $p = $Matrix[$i]->[$k]/$sum;
		$entropy1 -= $p*log($p)/$log2;
	}
	$info += $maxEntropy - $entropy1;
}
$info = int(1000*$info)/1000;
print OUT ">$pfm_name\t$pattern\t$patternRev\t$freq\t$ratio\t$info\t$score\t$FDR\t$palindrome\t$Nmembers\t$method\n";
my @items = ($pfm_name,$pattern,$patternRev,$freq,$ratio,$info,$score,$FDR,$palindrome,$Nmembers,$method);
print "<h3>Motif $pfm_name saved to $file_motif_save</h3>\n";
print "<TABLE BORDER=1><TR>\n";
for(my $j=0; $j<@headersGlob; ++$j){
	if($headerInd[$j]<0){ next; }
	print "<TD width=80><center>$headersGlob[$j]";
}
print "<TR>\n";
for(my $j=0; $j<@headersGlob; ++$j){
	if($headerInd[$j]<0){ next; }
	print "<TD><center>$items[$headerInd[$j]]\n";
}
print "</TABLE>\n";
print "<h3>Position frequency matrix (PFM)</h3>\n";
print "<TABLE BORDER=0><TR><TD><b><font color=blue>Position<TD width=70><center><b><font color=blue>A<TD width=70><center><b><font color=blue>C<TD width=70><center><b><font color=blue>G<TD width=70><center><b><font color=blue>T\n";
for(my $i=0; $i<$lenNew; ++$i){
	print OUT $i;
	print "<TR><TD>$i";
	for(my $k=0; $k<4; ++$k){
		my $x = $Matrix[$i]->[$k];
		print OUT "\t$x";
		print "<TD><center>$x\n";
	}
	print OUT "\n";
	print "\n";
}
print OUT "\n";
close OUT;
print "</TABLE>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  generate_PFM
#**************************************
{
my $outputFind = $hashInput{"file_motif_elem"};
my $file_motif_output = $hashInput{"file_motif_output"};
my $cluster = $hashInput{"cluster"};
my $pfm_start = $hashInput{"pfm_start"};
my $pfm_end = $hashInput{"pfm_end"};
my $pfm_name = $hashInput{"pfm_name"};
my $file_motif_save = $hashInput{"file_motif_save"};
my $reverse_PFM = $hashInput{"reverse_PFM"};

my $save_motif = 1;
if(!$file_motif_save){
	error_message("File motif is empty!");
}
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
my $found=0;
if(!open(OUT, ">$PATH_INFO/$loginname-config1.txt")){
	error_message("Cannot update configuration file!");
}
my $line;
while($line = <INFO>){
	my %hash=();
	if(length($line)<3){ next; }
	read_config_line($line,\%hash);
	my $file_motif1 = $hash{"type_motif"};
	if($file_motif1 eq $file_motif_save){
		$found=1;
	}
	print OUT $line;
}
if(!$found){
	print OUT "type_motif=$file_motif_save\tdescription=Motifs new\n";
}
close INFO;
close OUT;
if(!$found){
	system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
}
system("rm $PATH_INFO/$loginname-config1.txt");
if(!open(INFO, "$PATH_OUTPUT/$file_motif_output")){
	error_message("Cluster output file not found!");
}
if(!$found){
	open(OUT, ">$PATH_DATA/$loginname-$file_motif_save");
}else{
	open(INFO1, "$PATH_DATA/$loginname-$file_motif_save");
	my $last_line;
	while($line = <INFO1>){
		if($line =~ /^>$pfm_name\t/){
			error_message("Motif \'$pfm_name\' is present in file \'$file_motif_save\'! Use different name");
		}
		$last_line = $line;
	}
	close INFO1;
	open(OUT, ">>$PATH_DATA/$loginname-$file_motif_save");
	if(length($last_line) > 2){
		print OUT "\n";
	}
}
close INFO;

open(INFO, "$PATH_OUTPUT/$file_motif_output");
my $line;
my @headers;
my @headerInd;
my $junk;
while($line = <INFO>){
	if($line =~ /^Headers|^Parameters/){
		print OUT $line;
		if($line =~ /^Headers/){
			$line =~ s/\n$//;
			($junk,@headers) = split("\t",$line);
		}
	}	
	if($line =~ /^>$cluster\t/){ last; }
}
for(my $j=0; $j<@headersGlob; ++$j){
	$headerInd[$j]=-1;
	for(my $i=0; $i<@headers; ++$i){
		if($headersGlob[$j] eq $headers[$i]){ $headerInd[$j]=$i; }
	}
}
if($line !~ /^>$cluster\t/){ error_message("Cluster $cluster not found"); }
$line =~ s/^>//;
$line =~ s/\n$//;
my ($name,$pattern,$patternRev,$freq,$ratio,$info,$score,$FDR,$palindrome,$Nmembers,$method)=split(/\t/,$line);

# Print page header
print "<HTML><HEAD><TITLE>Motif $cluster</TITLE>\n";
print_header();
print "<h3>Motif \'$pfm_name\' is saved</h3>\n";
while($line = <INFO>){
	if($line =~ /^>/ || $line =~ /^members/){ last; }
}
if(!$line || $line =~ /^>/){  print "Member motifs not found"; exit(0); }
my %hashMotif;
while($line = <INFO>){
	$line =~ s/\n$//;
	if(!$line){ last; }
	my ($name,$pattern,$freq1,$info,$dir,$offset) = split("\t",$line);
	if(!$freq1){ $freq1=100; }
	if($hashInput{"use_motif_$name"}){
		$hashMotif{$name} = [$dir,$offset,$freq1];
	}
}
close INFO;
if(!%hashMotif){ print "Member motifs not found"; exit(0); }
if(!open(INFO, "$PATH_OUTPUT/$outputFind")){
	print "Cluster output file not found!\n"; exit(0);
}
my @sumFreq=();
my $nMotif=0;
my @Matrix=();
my $len = $pfm_end-$pfm_start+1;
for(my $i=0; $i<$len; ++$i){ $Matrix[$i]=[0,0,0,0]; }
while($line = <INFO>){
	$line =~ s/\n$//;
	if($line !~ /^>/){ next; }
	my ($name,$pattern) = split("\t",$line);
	$name =~ s/^>//;
	my $ref = $hashMotif{$name};
	if(!$ref){ next; }
	my($dir,$offset,$freq1) = @$ref;
	my ($name,@items) = split(/\t/,$line);
	my @Matrix1=();
	while($line = <INFO>){
		$line =~ s/\n$//;
		if(!$line){ last; }
		my ($pos,@m) = split(/\t/,$line);
		push(@Matrix1,\@m);
	}
	my $len1 = @Matrix1;
	for(my $i=0; $i<$len1; ++$i){
		my $pos1 = $i+$offset-$pfm_end+1;
		if($dir<0){ $pos1 = $offset+$len1-$i-1-$pfm_end+1; }
		if($pos1<0 || $pos1>=$len){ next; }
		$sumFreq[$pos1] += $freq1;
		for(my $k=0; $k<4; ++$k){
			$Matrix[$pos1]->[$k] += $freq1*$Matrix1[$i]->[$k];
		}
	}
}
close INFO;
if($reverse_PFM){
	@sumFreq = reverse(@sumFreq);
	@Matrix = reverse(@Matrix);
	for(my $i=0; $i<$len; ++$i){ @{$Matrix[$i]} = reverse(@{$Matrix[$i]}); }
}
while(!$sumFreq[0] && @sumFreq){
	shift(@sumFreq);
	shift(@Matrix);
}
$len = @Matrix;
$pattern="";
$info=0;
my $maxEntropy = 2;
my $log2 = log(2);
for(my $i=0; $i<$len; ++$i){
	if(!$sumFreq[$i]){ last; }
	my $max=0;
	my $sum=0;
	for(my $k=0; $k<4; ++$k){
		if($max < $Matrix[$i]->[$k]){ $max = $Matrix[$i]->[$k]; }
		$sum += $Matrix[$i]->[$k];
	}
	my $x=0;
	my $coef=1;
	my $entropy1=0;
	for(my $k=0; $k<4; ++$k){
		my $y = 1;
		if($Matrix[$i]->[$k] < 0.5*$max){ $y=0; }
		$x += $y*$coef;
		$coef *= 2;
		my $p = $Matrix[$i]->[$k]/$sum;
		$entropy1 -= $p*log($p)/$log2;
		$Matrix[$i]->[$k] = int(100*$p+0.5);
	}
	$info += $maxEntropy - $entropy1;
	$pattern .= $codes[$x];
}
$info = int(1000*$info)/1000;
$patternRev = reverse($pattern);
$patternRev =~ tr/ACMGRSVTWYHKDBN/TGKCYSBAWRDMHVN/;
$Nmembers = keys %hashMotif;
if($save_motif){
	print OUT ">$pfm_name\t$pattern\t$patternRev\t$freq\t$ratio\t$info\t$score\t$FDR\t$palindrome\t$Nmembers\t$method\n";
}
my @items = ($pfm_name,$pattern,$patternRev,$freq,$ratio,$info,$score,$FDR,$palindrome,$Nmembers,$method);
print "<TABLE BORDER=1><TR>\n";
for(my $j=0; $j<@headersGlob; ++$j){
	if($headerInd[$j]>=0){ print "<TD width=80><center>$headersGlob[$j]"; }
}
print "<TR>\n";
for(my $j=0; $j<@headersGlob; ++$j){
	if($headerInd[$j]<0){ next; }
	print "<TD><center>$items[$headerInd[$j]]\n";
}
print "</TABLE>\n";
print "<TABLE BORDER=1><TR><TD><b><font color=blue>Position<TD width=70><center><b><font color=blue>A<TD width=70><center><b><font color=blue>C<TD width=70><center><b><font color=blue>G<TD width=70><center><b><font color=blue>T\n";
for(my $i=0; $i<$len; ++$i){
	if($save_motif){ print OUT $i; }
	print "<TR><TD>$i";
	for(my $k=0; $k<4; ++$k){
		my $x = $Matrix[$i]->[$k];
		if($save_motif){ print OUT "\t$x"; }
		print "<TD><center>$x\n";
	}
	if($save_motif){ print OUT "\n"; }
	print "\n";
}
if($save_motif){
	print OUT "\n";
	close OUT;
}
print "</TABLE>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  motif_improve
#**************************************
{
my %hashConfig=%hashInput;
if(!$hashInput{"data_saved"}){
	if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){ error_message("Profile not found"); }
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash=();
		read_config_line($line,\%hash);
		if($hashInput{"mainPage"}){
			if($hash{"mainPage"}){ %hashConfig = %hash; last; }
		}else{
			if($dataname && $dataname eq $hash{"data"}){ %hashConfig = %hash; last; }
		}
	}
	close INFO;
	if($hashInput{"mainPage"}){ $dataname=""; }
}
my $file_motif = $hashConfig{"file_motif1"};
my $file_test = $hashConfig{"file_test"};
my $file_control = $hashConfig{"file_control"};
my $use_repeats = $hashConfig{"use_repeats"};
my $strand = $hashConfig{"strand"};
my $presence = $hashConfig{"presence"};

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<h3>Specify parameters of motif improvement</h3>\n";

print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"file_test\" TYPE=\"hidden\" VALUE=\"$file_test\">\n";
print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=\"$file_motif\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
if($hashInput{"data_saved"}){
	print "<INPUT NAME=\"data_saved\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
}
print "<h3>Files selected for motif improvement</h3>\n";
print "<TABLE border=1>\n";
print "<tr><td><b>File type<td><b>File name\n";
print "<tr><td><b>Motif file<td><font color=brown>$file_motif\n";
print "<tr><td><b>Sequence file #1 (test)<td><font color=brown>$file_test\n";
if($file_control){
	print "<tr><td><b>Sequence file #2 (control)<td><font color=brown>$file_control\n";
	print "<INPUT NAME=\"file_control\" TYPE=\"hidden\" VALUE=\"$file_control\">\n";
}
print "</TABLE>\n";

print "<h3>Parameters for motif improvement (modify if needed)</h3>\n";
my @yesno = ("No","Yes");
my @strand_list = ("Both","Forward","Forward vs. back");
my @iteration_list = (1,2,3,4,5);
my @subiteration_list = (1,3,5);
my @falsePositives = (0.1,0.2,0.5,1,2,3,5,7,10,15,20);
my @methods = ("none","Regression","Simple Resampling","Subtraction");
print "<TABLE BORDER=0>\n";
print "<TR><TD><select name=\"use_repeats\">\n";
for(my $i=0; $i<2; ++$i){ 
	print "<option value=$i";
	if($use_repeats==$i){ print " selected"; }
	print "> $yesno[$i]\n";
}
print "</select><TD>Use repeats for search\n";
print "<TR><TD><select name=\"strand\">\n";
for(my $i=0; $i<@strand_list; ++$i){ 
	print "<option value=$i";
	if($strand==$i){ print " selected"; }
	print ">$strand_list[$i]\n";
}
print "</select><TD>Search in strands\n";
print "<TR><TD><select name=\"falsePositives\">\n";
for(my $i=0; $i<@falsePositives; ++$i){ 
	print "<option value=$falsePositives[$i]";
	if($falsePositives[$i]==5){ print " selected"; }
	print ">$falsePositives[$i]\n";
}
print "</select><TD>Number of false positive hits per 10 Kb\n";
print "<TR><TD><select name=\"resampleMethod\">\n";
for(my $i=1; $i<@methods; ++$i){ 
	print "<option value=$i>$methods[$i]\n";
}
print "</select><TD>Resample method\n";
print "<TR><TD><select name=\"iterations\">\n";
for(my $i=0; $i<@iteration_list; ++$i){ 
	print "<option value=$iteration_list[$i]";
	if($iteration_list[$i]==3){ print " selected"; }
	print ">$iteration_list[$i]\n";
}
print "</select><TD>Select the number of iterations\n";
print "<TR><TD><select name=\"subiterations\">\n";
for(my $i=0; $i<@subiteration_list; ++$i){ 
	print "<option value=$subiteration_list[$i]";
	if($subiteration_list[$i]==5){ print " selected"; }
	print ">$subiteration_list[$i]\n";
}
print "</select><TD>Select the number of sub-iterations\n";
print "<TR><TD><select name=\"presence\">\n";
for(my $i=0; $i<2; ++$i){ 
	print "<option value=$i";
	if($presence==$i){ print " selected"; }
	print ">$yesno[$i]\n";
}
print "</select><TD>Count motif once per sequence\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"motif_improve1\">\n";
if($hashInput{"mainPage"}){
	print "<INPUT NAME=\"mainPage\" TYPE=\"hidden\" VALUE=\"mainPage\">\n";
}
print "<TR><TD><INPUT NAME=\"submit\" TYPE=submit VALUE=Submit>\n";
print "</FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  motif_improve1
#**************************************
{
my $file_motif = $hashInput{"file_motif"};
my $file_test = $hashInput{"file_test"};
my $file_control = $hashInput{"file_control"};
my $use_repeats = $hashInput{"use_repeats"};
my $strand = $hashInput{"strand"};
my $presence = $hashInput{"presence"};
my $iterations = $hashInput{"iterations"};
my $subiterations = $hashInput{"subiterations"};
my $resampleMethod = $hashInput{"resampleMethod"};
my $falsePositives = $hashInput{"falsePositives"};

my $fileMotif = $file_motif;
if($file_motif !~ /^public-/ && !$hashInput{"data_saved"}){ $fileMotif = "$loginname-$file_motif"; }
my $fileTest = $file_test;
if($file_test !~ /^public-/ && !$hashInput{"data_saved"}){ $fileTest = "$loginname-$file_test"; }
my $fileControl = $file_control;
if($file_control && $file_control !~ /^public-/ && !$hashInput{"data_saved"}){ $fileControl = "$loginname-$file_control"; }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";

my @motif_list=();
if(open(INFO, "$PATH_INFO/$loginname-config.txt")){
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash1=();
		read_config_line($line,\%hash1);
		if($hash1{"type_motif"}){ push(@motif_list,$hash1{"type_motif"}); }
	}
	close INFO;
	@motif_list = sort @motif_list;
}
my $list1 = join("\",\"",@motif_list);
if($list1){ $list1 = "\"".$list1."\""; }
print "motif_list = new Array($list1);\n";
print "function save_motif_file(){\n";
print "	if(!document.cisfinder.file_motif.value){\n";
print "		alert(\"Enter name for new motif file\"); return(false);\n";
print "	}\n";
print "	var filename = document.cisfinder.file_motif.value;\n";
print "	if(filename.search(/ |\\/|\\&|\\t|\\%|\\+|\\@|\\!/) >= 0){\n";
print "		alert(\"File name should have neither spaces nor special characters\");\n";
print "		return(false);\n";
print "	}\n";
print "	var i;\n";
print "	for(i=0; i<motif_list.length; ++i){\n";
print "		if(filename == motif_list[i] && !confirm(\"Do you want to overwrite existing file \"+filename+\"?\")){\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	document.cisfinder.action.value=\"save_motifs\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function show_motifs(){\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value=\"show_motifs\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();

# Assemble the command line and execute patternTest
my $outputIDtest = get_outputID(2);
my $outputIDprog = $outputIDtest+1;
my @command=(
	"./patternTest",
	"-i $PATH_DATA/$fileMotif",
	"-f $PATH_DATA/$fileTest",
	"-o $PATH_OUTPUT/$outputIDtest.txt",
	"-prog $PATH_OUTPUT/$outputIDprog.txt",
	"-iter $iterations",
	"-siter $subiterations",
	"-method $resampleMethod",
	"-fp $falsePositives",
);
if($file_control){ push(@command,"-c $PATH_DATA/$fileControl"); }
if($use_repeats){  push(@command,"-userep"); }
if($strand){       push(@command,"-strand $strand"); }
if($presence){     push(@command,"-one"); }

print "<b>Improving motifs...</b><br>\n";
my $string = join(' ',@command);
#print "$string\n";
my $response = `$string`;
print "<pre>$response</pre>\n";
if($response =~ /ERROR/i){
	print "<h3>Errors in program run<h3>Output was not generated<p>\n";
	exit(0);
}
if(!open(INFO, "$PATH_OUTPUT/$outputIDtest.txt")){
	print "Motif output file not found!\n"; exit(0);
}
my $count=0;
while(my $line = <INFO>){
	if($line =~ /^>/){ ++$count; }
}
close INFO;
print "Number of motifs = $count<br>\n";
print "Output as plain text: <a href=../output/$outputIDtest.txt target=_BLANK483>Motif file</a><br>\n";
print "View progress: <a href=../output/$outputIDprog.txt target=_BLANK493>Progress file</a><p>\n";
print "<INPUT NAME=\"show_motif_button\" TYPE=button VALUE=\"Show improved motifs\" onClick=show_motifs();><p>\n";

# Save new motif file
$file_motif =~ s/public-//;
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST\">\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"save_motifs\">\n";
print "<INPUT NAME=\"file_motif_output\" TYPE=\"hidden\" VALUE=$outputIDtest.txt>\n";
print "<b>Save the file with improved motifs</b><br>\n";
print "To overwrite the original file, specify the same name<br>\n";
print "File name: <INPUT NAME=\"file_motif\" SIZE=15 VALUE=$file_motif-improved>\n";
print "Description: <INPUT NAME=\"description_motif\" SIZE=50 VALUE=\"Motifs from $file_motif over-represented in sequences $file_test\"> (modify if needed)\n";
print "<INPUT NAME=\"save_motif_button\" TYPE=button VALUE=\"Save motifs\" onClick=\"return save_motif_file();\">\n";
print "<INPUT NAME=\"close_button\" TYPE=button VALUE=\"Cancel\" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  motif_compare
#**************************************
{
my %hashConfig=%hashInput;
my @motif_list;
my $file_motif1;
if(!$hashInput{"data_saved"}){
	if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){ error_message("Profile not found"); }
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash=();
		read_config_line($line,\%hash);
		my $file_motif = $hash{"type_motif"};
		if($file_motif){ push(@motif_list,[$file_motif,$hash{"description"}]); next; }
		if($hashInput{"mainPage"}){
			if($hash{"mainPage"}){ %hashConfig = %hash; }
		}else{
			if($dataname && $dataname eq $hash{"data"}){ %hashConfig = %hash; }
		}
	}
	close INFO;
	if($hashInput{"mainPage"}){ $dataname=""; }
	@motif_list = sort {lc($a->[0]) cmp lc($b->[0])} @motif_list;
	$file_motif1 = $hashConfig{"file_motif1"};
}
if($loginname ne "public" && open(INFO, "$PATH_INFO/public-config.txt")){
	my @motif_list1=();
	while(my $line = <INFO>){
		my %hash=();
		read_config_line($line,\%hash);
		my $file_motif = $hash{"type_motif"};
		if($file_motif){ push(@motif_list1,["public-$file_motif",$hash{"description"}]); }
	}
	push(@motif_list, sort {lc($a->[0]) cmp lc($b->[0])} @motif_list1);
}
close INFO;
my $file_motif1 = $hashConfig{"file_motif1"};
my $match_thresh = $hashConfig{"match_thresh"};
if(!$match_thresh){ $match_thresh=0.75; }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
my ($items,$descriptions) = get_array_lists(\@motif_list);
print "motif_list = new Array($items);\n";
print "motif_description = new Array($descriptions);\n";
print "function update_description() {\n";
print "  index = document.cisfinder.file_motif2.selectedIndex;\n";
print "  if(index>=0){ document.cisfinder.description_motif2.value = motif_description[index]; }\n";
print "  else{ document.cisfinder.description_motif2.value = \"\"; }\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header("update_description();");
print "<h3>Comparing motifs (file $file_motif1)</h3>\n";

print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST >\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"motif_compare1\">\n";
print "<INPUT NAME=\"file_motif1\" TYPE=\"hidden\" VALUE=\"$file_motif1\">\n";
#if($hashInput{"data_saved"}){
#	print "<INPUT NAME=\"data_saved\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
#}
if($hashInput{"mainPage"}){
	print "<INPUT NAME=\"mainPage\" TYPE=\"hidden\" VALUE=\"mainPage\">\n";
}
elsif($dataname){
	print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
}

print "<b>Select motif file for comparison:</b><br>\n";
print "<TABLE border=0>\n";
print "<TR><TD><center><b>File<TD><center><b>Description\n";
print "<TR><TD><select name=\"file_motif2\" onChange=update_description();>\n";
for(my $i=0; $i<@motif_list; ++$i){ 
	print "<option value=\"$motif_list[$i]->[0]\"";
	if($motif_list[$i]->[0] =~ /CisView/){ print " selected"; }
	print "> $motif_list[$i]->[0]\n";
}
print "</select><TD><INPUT NAME=\"description_motif2\" SIZE=40>\n";
print "</TABLE><p>\n";
print "<b>Select match threshold (correlation): </b>\n";
my @match_thresh = (0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.98);
print "<TR><TD><select name=\"match_thresh\">\n";
for(my $i=0; $i<@match_thresh; ++$i){ 
	print "<option value=$match_thresh[$i]";
	if($match_thresh[$i]==0.75){ print " selected"; }
	print "> $match_thresh[$i]\n";
}
print "</select><p>\n";
print "</FORM><p>\n";
print "<INPUT NAME=\"submit_button\" TYPE=button VALUE=\"Compare motifs\" LANGUAGE=\"javascript\" onClick=\"document.cisfinder.submit();\"><p>\n";
print "<INPUT NAME=\"close_button\" TYPE=button VALUE=\"Cancel\" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  motif_compare1
#**************************************
{
if($hashInput{"mainPage"}){ $dataname=""; }
my $file_motif1 = $hashInput{"file_motif1"};
my $file_motif2 = $hashInput{"file_motif2"};
my $match_thresh = $hashInput{"match_thresh"};
if(!$match_thresh || !$file_motif1 || !$file_motif2 ){ error_message("ERR in motif_compare!"); }
my $fileMotif1 = $file_motif1;
if($fileMotif1 !~ /^public-/ && !$hashInput{"data_saved"}){ $fileMotif1 = "$loginname-$fileMotif1"; }
my $fileMotif2 = $file_motif2;
if($fileMotif2 !~ /^public-/ && !$hashInput{"data_saved"}){ $fileMotif2 = "$loginname-$fileMotif2"; }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "function show_motif(file,motif){\n";
print "	document.cisfinder.motif.value = motif;\n";
print "	document.cisfinder.file_motif.value = file;\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<h3>Comparing motif files $file_motif1 and $file_motif2</h3>\n";

print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST >\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"show_motif\">\n";
print "<INPUT NAME=\"motif\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=\"\"></FORM>\n";

# Assemble the command line and execute
my $outputIDcompare = get_outputID(1);
my @command=(
	"./patternCompare",
	"-i1 $PATH_DATA/$fileMotif1",
	"-i2 $PATH_DATA/$fileMotif2",
	"-o $PATH_OUTPUT/$outputIDcompare.txt",
	"-match $match_thresh"
);
print "<b>Comparing motifs...</b><br>\n";
my $string = join(' ',@command);
#print "$string\n";
my $response = `$string`;
print "<pre>$response</pre>\n";
if($response =~ /ERROR/i){
	print "<h3>Errors in program run<h3>Output was not generated<p>\n";
	exit(0);
}
if(!open(INFO, "$PATH_OUTPUT/$outputIDcompare.txt")){
	print "Motif output file not found!\n"; exit(0);
}
my $count=0;
my @lines=();
while(my $line = <INFO>){
	if($line =~ /^>/){ ++$count; }
	$line =~ s/\n$//;
	push(@lines,$line);
}
close INFO;
print "Number of motifs compared = $count<br>\n";
print "Output as plain text: <a href=../output/$outputIDcompare.txt target=_BLANK642>Motif compare file</a><p>\n";
print "<HR NOSHADE COLOR=BLUE></HR>\n";

print "<h3>Table of aligned motifs</h3>\n";
print "<TABLE BORDER=0><TR><TD>No.<TD><b><font color=blue>Motif Name<TD><center><b><font color=blue>Pattern<TD><b><font color=blue>Match(r)<TD><b><font color=blue>Strand<TD><b><font color=blue>Offset\n";
my $count=1;
my $countMatch = 0;
foreach my $line (@lines){
	$line =~ s/\n$//;
	if(!$line){ next; }
	my ($name,$pattern,$match,$dir,$offset) = split("\t",$line);
	if($name =~ /^>/){
		$name =~ s/^>//;
		print "<TR><TD>$count<TD bgcolor=\"#FFFF80\"><b>";
		print "<a href=\"JavaScript: show_motif('$file_motif1','$name');\">$name</a>";
		print "<TD><pre><font color=brown>$pattern\n";
		++$count;
		$countMatch=0;
	}else{
		++$countMatch;
		if($countMatch>=10){ next; }
		print "<TR><TD><TD><font size=-1>";
		print "<a href=\"JavaScript: show_motif('$file_motif2','$name');\">$name</a>";
		print "<TD><pre>$pattern<TD><center>$match<TD><center>$dir<TD><center>$offset\n";
	}
}
print "</TABLE>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  motif_cluster
#**************************************
{
my %hashConfig=%hashInput;
if(open(INFO, "$PATH_INFO/$loginname-config.txt")){
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash=();
		read_config_line($line,\%hash);
		if($hashInput{"mainPage"}){
			if($hash{"mainPage"}){ %hashConfig = %hash; }
		}else{
			if($dataname && $dataname eq $hash{"data"}){ %hashConfig = %hash; }
		}
	}
	close INFO;
	if($hashInput{"mainPage"}){ $dataname=""; }
}
my $file_motif1 = $hashConfig{"file_motif1"};
my $match_thresh = $hashConfig{"match_thresh"};
if(!$match_thresh){ $match_thresh = 0.75; }
if(!$file_motif1){ error_message("ERR in motif_cluster!"); }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<h3>Clustering motifs from file $file_motif1</h3>\n";

print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST >\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"motif_cluster1\">\n";
print "<INPUT NAME=\"file_motif1\" TYPE=\"hidden\" VALUE=\"$file_motif1\">\n";
if($hashInput{"data_saved"}){
	print "<INPUT NAME=\"data_saved\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
}
if($hashInput{"mainPage"}){
	print "<INPUT NAME=\"mainPage\" TYPE=\"hidden\" VALUE=\"mainPage\">\n";
}
elsif($dataname){
	print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
}

print "<b>Select match threshold (correlation) for clustering:</b><br>\n";
my @match_thresh = (0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.98);
print "<TR><TD><select name=\"match_thresh\">\n";
for(my $i=0; $i<@match_thresh; ++$i){ 
	print "<option value=$match_thresh[$i]";
	if($match_thresh[$i]==0.75){ print " selected"; }
	print "> $match_thresh[$i]\n";
}
print "</select><p>\n";
print "</FORM><p>\n";
print "<INPUT NAME=\"submit_button\" TYPE=button VALUE=\"Cluster motifs\" LANGUAGE=\"javascript\" onClick=\"document.cisfinder.submit();\"><p>\n";
print "<INPUT NAME=\"close_button\" TYPE=button VALUE=\"Cancel\" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  motif_cluster1
#**************************************
{
if($hashInput{"mainPage"}){ $dataname=""; }
my $file_motif1 = $hashInput{"file_motif1"};
my $match_thresh = $hashInput{"match_thresh"};
if(!$match_thresh || !$file_motif1){ error_message("ERR in motif_cluster1!"); }
my $fileMotif1 = $file_motif1;
if($fileMotif1 !~ /^public-/ && !$hashInput{"data_saved"}){ $fileMotif1 = "$loginname-$fileMotif1"; }
my @motif_list=();
if(open(INFO, "$PATH_INFO/$loginname-config.txt")){
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash1=();
		read_config_line($line,\%hash1);
		if($hash1{"type_motif"}){ push(@motif_list,$hash1{"type_motif"}); }
	}
	close INFO;
	@motif_list = sort @motif_list;
}
my $list1 = join("\",\"",@motif_list);
if($list1){ $list1 = "\"".$list1."\""; }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "motif_list = new Array($list1);\n";
print "function save_motif_file(){\n";
print "	if(!document.cisfinder.file_motif.value){\n";
print "		alert(\"Enter name for new motif file\"); return(false);\n";
print "	}\n";
print "	var filename = document.cisfinder.file_motif.value;\n";
print "	if(filename.search(/ |\\/|\\&|\\t|\\%|\\+|\\@|\\!/) >= 0){\n";
print "		alert(\"File name should have neither spaces nor special characters\");\n";
print "		return(false);\n";
print "	}\n";
print "	var i;\n";
print "	for(i=0; i<motif_list.length; ++i){\n";
print "		if(filename == motif_list[i] && !confirm(\"Do you want to overwrite existing file \"+filename+\"?\")){\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value=\"save_motifs\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function show_motifs(){\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value=\"show_motifs\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<h3>Clustering motifs from file $file_motif1</h3>\n";

# Assemble the command line and execute
my $outputIDclust = get_outputID(2);
my $outputIDfind = $outputIDclust+1;
my @command=(
	"./patternCluster",
	"-i $PATH_DATA/$fileMotif1",
	"-o $PATH_OUTPUT/$outputIDclust.txt",
	"-match $match_thresh"
);
print "<b>Clustering motifs...</b><br>\n";
my $string = join(' ',@command);
#print "$string\n";
my $response = `$string`;
print "<pre>$response</pre>\n";
if($response =~ /ERROR/i){
	print "<h3>Errors in program run<h3>Output was not generated<p>\n";
	exit(0);
}
`cp $PATH_DATA/$fileMotif1 $PATH_OUTPUT/$outputIDfind.txt`;
if(!open(INFO, "$PATH_OUTPUT/$outputIDclust.txt")){
	print "Cluster output file not found!\n"; exit(0);
}
my $count=0;
while(my $line = <INFO>){
	if($line =~ /^>/){ ++$count; }
}
close INFO;
my $nPages = int(($count+9)/10);
print "Number of clusters = $count<p>\n";

# Display results
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"show_motifs\">\n";
print "<INPUT NAME=\"page\" TYPE=\"hidden\" VALUE=1>\n";
print "<INPUT NAME=\"nPages\" TYPE=\"hidden\" VALUE=$nPages>\n";
print "<INPUT NAME=\"nClusters\" TYPE=\"hidden\" VALUE=$count>\n";
print "<INPUT NAME=\"file_motif_elem\" TYPE=\"hidden\" VALUE=$outputIDfind.txt>\n";
print "<INPUT NAME=\"file_motif_output\" TYPE=\"hidden\" VALUE=$outputIDclust.txt>\n";
print "<INPUT NAME=\"show_cluster_button\" TYPE=button VALUE=\"Show clusters of motifs\" onClick=show_motifs();>\n";
if($hashInput{"data_saved"}){
	print "<INPUT NAME=\"data_saved\" TYPE=\"hidden\" VALUE=1>\n";
}
my $file_motif_save = "$file_motif1-clusters";
$file_motif_save =~ s/^$loginname-|^public-//;
print "File name: <INPUT NAME=\"file_motif\" SIZE=15 VALUE=$file_motif_save>\n";
print "Description: <INPUT NAME=\"description_motif\" SIZE=50 VALUE=\"$file_motif1-clusters\"> (modify if needed)<br>\n";
print "<INPUT NAME=\"save_motif_button\" TYPE=button VALUE=\"Save motifs\" onClick=save_motif_file();>\n";

print "</FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  motif_search
#**************************************
{
my %hashConfig=%hashInput;
if(!$hashInput{"data_saved"}){
	if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){ error_message("Profile not found"); }
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash=();
		read_config_line($line,\%hash);
		if($hashInput{"mainPage"}){
			if($hash{"mainPage"}){ %hashConfig = %hash; last; }
		}else{
			if($dataname && $dataname eq $hash{"data"}){ %hashConfig = %hash; last; }
		}
	}
	close INFO;
	if($hashInput{"mainPage"}){ $dataname=""; }
}
my $file_motif = $hashConfig{"file_motif1"};
my $file_test = $hashConfig{"file_test"};
my $use_repeats_search = $hashConfig{"use_repeats_search"};
my $use_redundant = $hashConfig{"use_redundant"};
my $strand = $hashConfig{"strand"};
my $presence = $hashConfig{"presence"};
my $interval = $hashConfig{"interval"};
if(!$interval){ $interval=100; }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder</TITLE>\n";
print_header();
print "<H3>Parameters of motif search</H3>\n";

my @yesno = ("No","Yes");
my @strand_list = ("Both","Forward");
my @falsepos_list = (1,2,3,5,7,10,15,20);
my @addThresh_list = (-3,-2,-1,0,1,2,3);
my @interval_list = (5,10,20,50,100,200,500,1000,2000);
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD><font color=brown>$file_motif<TD>Motif file selected\n";
print "<TR><TD><font color=brown>$file_test<TD>Sequence file selected\n";
print "<TR><TD><select name=\"falsePositives\">\n";
print "<option value=0>Use existing thresholds";
for(my $i=0; $i<@falsepos_list; ++$i){ 
	print "<option value=$falsepos_list[$i]";
	if($falsepos_list[$i]==5){ print " selected"; }
	print ">$falsepos_list[$i]\n";
}
print "</select><TD>Number of false positives per 10 Kb\n";
print "<TR><TD><select name=\"addThreshold\">\n";
for(my $i=0; $i<@addThresh_list; ++$i){ 
	print "<option value=\"$addThresh_list[$i]\"";
	if(!$addThresh_list[$i]){ print " selected"; }
	print "> $addThresh_list[$i]\n";
}
print "</select><TD>Add this value to the score threshold\n";
print "<TR><TD><select name=\"strand\">\n";
for(my $i=0; $i<@strand_list; ++$i){ 
	print "<option value=$i>$strand_list[$i]\n";
}
print "</select><TD>Search in strands\n";
print "<TR><TD><select name=\"use_repeats_search\">\n";
print "<option value=0>No\n";
print "<option value=1";
if($use_repeats_search){ print " selected"; }
print ">Yes\n";
print "</select><TD>Search repeats\n";
print "<TR><TD><select name=\"use_redundant\">\n";
print "<option value=0>No\n";
print "<option value=1";
if($use_redundant){ print " selected"; }
print ">Yes\n";
print "</select><TD>Include overlapping and redundant motifs in the output\n";
print "<TR><TD><select name=\"interval\">\n";
for(my $i=0; $i<@interval_list; ++$i){ 
	print "<option value=$interval_list[$i]";
	if($interval==$interval_list[$i]){ print " selected"; }
	print "> $interval_list[$i]\n";
}
print "</select><TD>Interval (bp) for frequency distribution\n";

print "<TR><TD><INPUT NAME=\"submit\" TYPE=submit VALUE=Submit>\n";
print "</TABLE>\n";
my $taskID = get_outputID(4);
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"motif_search1\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=\"$file_motif\">\n";
print "<INPUT NAME=\"file_test\" TYPE=\"hidden\" VALUE=\"$file_test\">\n";
print "<INPUT NAME=\"taskID\" TYPE=\"hidden\" VALUE=\"$taskID\">\n";
if($hashInput{"data_saved"}){
	print "<INPUT NAME=\"data_saved\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
}
if($hashInput{"mainPage"}){
	print "<INPUT NAME=\"mainPage\" TYPE=\"hidden\" VALUE=\"mainPage\">\n";
}
print "</FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  motif_search1
#**************************************
{
my $file_motif = $hashInput{"file_motif"};
my $file_fasta = $hashInput{"file_test"};
my $falsePositives = $hashInput{"falsePositives"};
my $addThresh = $hashInput{"addThresh"};
my $strand = $hashInput{"strand"};
my $use_repeats_search = $hashInput{"use_repeats_search"};
my $use_redundant = $hashInput{"use_redundant"};
my $interval = $hashInput{"interval"};
my $taskID = $hashInput{"taskID"};

my $filenameMotif = "$loginname-$file_motif";
if($file_motif =~ /^public-/ || $hashInput{"data_saved"}){ $filenameMotif = $file_motif; }
my $fileFasta = "$loginname-$file_fasta";
if($file_fasta =~ /^public-/ || $hashInput{"data_saved"}){ $fileFasta = $file_fasta; }

my @search_list=();
if(open(INFO, "$PATH_INFO/$loginname-config.txt")){
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash1=();
		read_config_line($line,\%hash1);
		if($hash1{"type_search"}){ push(@search_list,$hash1{"type_search"}); }
	}
	close INFO;
	@search_list = sort @search_list;
}
my $list1 = join("\",\"",@search_list);
if($list1){ $list1 = "\"".$list1."\""; }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "search_list = new Array($list1);\n";
print "function save_search_file(){\n";
print "	if(!document.cisfinder.file_search.value){\n";
print "		alert(\"Enter name for new search results file\"); return(false);\n";
print "	}\n";
print "	var filename = document.cisfinder.file_search.value;\n";
print "	if(filename.search(/ |\\/|\\&|\\t|\\%|\\+|\\@|\\!/) >= 0){\n";
print "		alert(\"File name should have neither spaces nor special characters\");\n";
print "		return(false);\n";
print "	}\n";
print "	var i;\n";
print "	for(i=0; i<search_list.length; ++i){\n";
print "		if(filename == search_list[i] && !confirm(\"Do you want to overwrite existing file \"+filename+\"?\")){\n";
print "			return false;\n";
print "		}\n";
print "	}\n";
print "	document.cisfinder.search.value = 0;\n";
print "	document.cisfinder.action.value = \"save_search\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function show_motifs(){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"show_search_motifs\";\n";
print "	document.cisfinder.search.value = 0;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function show_sequences(){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"show_search_sequences\";\n";
print "	document.cisfinder.search.value = 0;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function search_sequences(){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"show_search_sequences\";\n";
print "	document.cisfinder.search.value = 1;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function get_info(iTF){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	var motif, ind;\n";
print "	if(iTF==1){ ind=document.cisfinder.TF1.selectedIndex; motif=document.cisfinder.TF1.options[ind].value; }\n";
print "	else if(iTF==2){ ind=document.cisfinder.TF2.selectedIndex; motif=document.cisfinder.TF2.options[ind].value; }\n";
print "	else if(iTF==3){ ind=document.cisfinder.TF3.selectedIndex; motif=document.cisfinder.TF3.options[ind].value; }\n";
print "	if(ind==0){\n";
print "		alert(\"Motif is not selected\");\n";
print "		return(false);\n";
print "	}\n";
print "	document.cisfinder.motif.value = motif;\n";
print "	document.cisfinder.action.value = \"show_motif\";\n";
print "	document.cisfinder.search.value = 1;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<b>Searching motifs ($file_motif) in sequences ($file_fasta)...</b><br>\n";
print "If connection is lost, use Refresh button<br>\n";

# Assemble the command line and execute patternScan
my @command=("./get_search.pl $loginname $sessionID $file_fasta $file_motif $fileFasta $filenameMotif $taskID");
if($use_repeats_search==1){ push(@command,"-userep"); }
if($use_redundant==1){ push(@command,"-redund"); }
if($addThresh){ push(@command,"-thresh $addThresh"); }
if($falsePositives){ push(@command,"-fp $falsePositives"); }
if($strand){ push(@command,"-strand $strand"); }
if($interval){ push(@command,"-int $interval"); }
if($dataname){ push(@command,"-data $dataname"); }
my $string = join(' ',@command);

my $response="";
if(open (INFO, "$PATH_OUTPUT/$taskID.txt")){
	while(my $line=<INFO>){ $response .= $line; }
	close INFO;
}else{
	#print "$string<p>\n";
	`$string`;
	$response="";
	if(open (INFO, "$PATH_OUTPUT/$taskID.txt")){
		while(my $line=<INFO>){ $response .= $line; }
		close INFO;
	}
}
if($response =~ /\<\/HTML\>/i){
	print $response;
}else{
	print "<b>Not finished yet... Try to refresh the page later.</b>\n";
}
exit(0);
}

#**************************************
sub  save_motifs
#**************************************
{
my $file_motif = $hashInput{"file_motif"};
my $file_motif_output = $hashInput{"file_motif_output"};
my $description_motif = $hashInput{"description_motif"};

open(OUT, ">$PATH_INFO/$loginname-config1.txt");
open(INFO, "$PATH_INFO/$loginname-config.txt");
while(my $line = <INFO>){
	if($line !~ /type_motif=$file_motif\t/){ print OUT $line; }
}
close INFO;
print OUT "type_motif=$file_motif\tdescription=$description_motif\n";
close OUT;
system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
system("rm $PATH_INFO/$loginname-config1.txt");
system("cp $PATH_OUTPUT/$file_motif_output $PATH_DATA/$loginname-$file_motif");

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print_header();
print "<h3>Motif file $file_motif saved</h3>\n";
print "<INPUT NAME=\"close_button\" TYPE=button VALUE=\"Close window\" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
return;
}

#**************************************
sub  save_search
#**************************************
{
my $file_search = $hashInput{"file_search"};
my $file_search_output = $hashInput{"file_search_output"};
my $description = $hashInput{"description"};

open(OUT, ">$PATH_INFO/$loginname-config1.txt");
open(INFO, "$PATH_INFO/$loginname-config.txt");
while(my $line = <INFO>){
	if($line !~ /type_search=$file_search\t/){ print OUT $line; }
}
close INFO;
print OUT "type_search=$file_search\tdescription=$description\n";
close OUT;
system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
system("rm $PATH_INFO/$loginname-config1.txt");
system("cp $PATH_OUTPUT/$file_search_output $PATH_DATA/$loginname-$file_search");

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print_header();
print "<h3>Search file $file_search saved</h3>\n";
print "You may need to push the 'Refresh' button in the main screen to see this file in the menu<p>\n";
print "<INPUT NAME=\"close_button\" TYPE=button VALUE=\"Close window\" LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  motif_delete
#**************************************
{
my $file_motif = $hashInput{"file_motif"};
my $file_motif_output = $hashInput{"file_motif_output"};
my $motif_delete = $hashInput{"motif_delete"};

my $fileMotif = "$PATH_DATA/$loginname-$file_motif";
if($file_motif_output){
	$fileMotif = "$PATH_OUTPUT/$file_motif_output";
}
if(!open(INFO, $fileMotif)){
	error_message("Motif file not found!");
}
if(!open(OUT, ">$fileMotif.bak")){
	error_message("Cannot update motif file!");
}
my $found = 0;
while(my $line = <INFO>){
	if($line =~ /^>$motif_delete/){
		$found = 1;
		while($line = <INFO>){
			if($line =~ /^>/){ last; }
		}
	}
	print OUT "$line";
}
close INFO;
close OUT;
if($found){
	system("cp $fileMotif.bak $fileMotif");
}
system("rm $fileMotif.bak");
return;
}

#**************************************
sub  motif_reverse
#**************************************
{
my $file_motif = $hashInput{"file_motif"};
my $file_motif_output = $hashInput{"file_motif_output"};
my $motif_reverse = $hashInput{"motif_reverse"};

my $fileMotif = "$PATH_DATA/$loginname-$file_motif";
if($file_motif_output){
	$fileMotif = "$PATH_OUTPUT/$file_motif_output";
}
if(!open(INFO, $fileMotif)){
	error_message("Motif file not found!");
}
if(!open(OUT, ">$fileMotif.bak")){
	error_message("Cannot update motif file!");
}
my $found = 0;
while(my $line = <INFO>){
	if($line =~ /^>$motif_reverse/){
		$found = 1;
		my ($name,$pattern,$patternRev,@items) = split(/\t/,$line);
		print OUT "$name\t$patternRev\t$pattern\t".join("\t",@items);
		my @Matrix=();
		my $count=0;
		while($line = <INFO>){
			$line =~ s/\n$//;
			if(!$line || $line =~ /^>/){ last; }
			my ($pos,@m) = split(/\t/,$line);
			if($pos != $count){ print "Error in matrix format!\n"; exit(0); }
			push(@Matrix,\@m);
			++$count;
		}
		@Matrix = reverse(@Matrix);
		for(my $i=0; $i<$count; ++$i){
			@{$Matrix[$i]} = reverse(@{$Matrix[$i]});
			print OUT $i;
			for(my $j=0; $j<4; ++$j){ print OUT "\t$Matrix[$i]->[$j]"; }
			print OUT "\n";
		}
		print OUT "\n";
		next;
	}
	print OUT "$line";
}
close INFO;
close OUT;
if($found){
	system("cp $fileMotif.bak $fileMotif");
}
system("rm $fileMotif.bak");
return;
}

#**************************************
sub  show_motifs
#**************************************
{
my $page = $hashInput{"page"};
my $nPages = $hashInput{"nPages"};
my $nMotifs = $hashInput{"nMotifs"};
my $file_motif_output = $hashInput{"file_motif_output"};
my $file_motif_elem = $hashInput{"file_motif_output"};
my $file_motif = $hashInput{"file_motif1"};
if(!$file_motif){ $file_motif = $hashInput{"file_motif"}; }

my $filenameMotif = "$PATH_DATA/$loginname-$file_motif";
if($file_motif =~ /^public-/ || $hashInput{"data_saved"}){ $filenameMotif = "$PATH_DATA/$file_motif"; }
if($file_motif_output){
	$filenameMotif = "$PATH_OUTPUT/$file_motif_output";
}elsif(!$file_motif){
	error_message("No motif file!");
}
if(!$page || !$nMotifs){
	$page=1;
	if(!open(INFO, $filenameMotif)){
		error_message("Motif file $filenameMotif not found!");
	}
	my $count=0;
	while(my $line = <INFO>){
		if($line =~ /^>/){ ++$count; }
	}
	close INFO;
	$nMotifs = $count;
	$nPages = int(($count+9)/10);
}
my $outputIDmotif = get_outputID(1);
system("cp $filenameMotif $PATH_OUTPUT/$outputIDmotif.txt");

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "var Npages = $nPages;\n";
print "function cleanAll(){\n";
print "	document.cisfinder.motif_delete.value = \"\";\n";
print "	document.cisfinder.motif_reverse.value = \"\";\n";
print "	document.cisfinder.target = \"\";\n";
print "	document.cisfinder.action.value = \"show_motifs\";\n";
print "}\n";
print "function nextpage(page){\n";
print "	cleanAll();\n";
print "	document.cisfinder.page.value = page;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function jump_to_page(){\n";
print "	var page = document.cisfinder.go_to_page.value;\n";
print "	if(isNaN(page)){ alert(\"Wrong page number!\"); return false; }\n";
print "	page = Math.ceil(page);\n";
print "	if(page<1 || page>Npages){ alert(\"Wrong page number!\"); return false; }\n";
print "	nextpage(page);\n";
print "}\n";
print "function search_motif(){\n";
print "	if(!document.cisfinder.search_term.value){ alert(\"Enter search term!\"); return false; }\n";
print "	cleanAll();\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"search_motif_name\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function show_motif(motif,logo){\n";
print "	cleanAll();\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.motif.value = motif;\n";
print "	document.cisfinder.seqLogo.value = logo;\n";
if($file_motif_elem){
	print "	document.cisfinder.action.value = \"show_cluster\";\n";
}else{
	print "	document.cisfinder.action.value = \"show_motif\";\n";
}
print "	document.cisfinder.submit();\n";
print "}\n";
print "function motif_delete(motif){\n";
print "	cleanAll();\n";
print "	document.cisfinder.motif_delete.value = motif;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function motif_reverse(motif){\n";
print "	cleanAll();\n";
print "	document.cisfinder.motif_reverse.value = motif;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();

# Generate motif-logo
my $start = ($page-1)*10+1;
my $outputIDlogo = get_outputID(10);
my $ID = $outputIDlogo;
for(my $imotif=$start; $imotif<$start+10 && $imotif<=$nMotifs; ++$imotif){
	`./motiflogo $filenameMotif $PATH_OUTPUT/$ID.gif -num $imotif`;
	++$ID;
}
print "<h3>Motifs in file <font color=BROWN>$file_motif</font> (P. $page out of $nPages)</h3>\n";
print "<b>Output in plain text:</b>\n";
print "<a href=../output/$outputIDmotif.txt target=_BLANK151>output</a><p>\n";
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<INPUT NAME=\"page\" TYPE=\"hidden\" VALUE=$page>\n";
print "<INPUT NAME=\"nPages\" TYPE=\"hidden\" VALUE=$nPages>\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"motif\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"nMotifs\" TYPE=\"hidden\" VALUE=\"$nMotifs\">\n";
print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=$file_motif>\n";
print "<INPUT NAME=\"motif_delete\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"motif_reverse\" TYPE=\"hidden\" VALUE=\"\">\n";
if($file_motif_output){
	print "<INPUT NAME=\"file_motif_output\" TYPE=\"hidden\" VALUE=\"$file_motif_output\">\n";
}
if($file_motif){
	print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=\"$file_motif\">\n";
}
if($file_motif_elem){
	print "<INPUT NAME=\"file_motif_elem\" TYPE=\"hidden\" VALUE=\"$file_motif_elem\">\n";
}
print "<INPUT NAME=\"seqLogo\" TYPE=\"hidden\" VALUE=\"\">\n";
my $newpage = $page-1;
if($page > 1){ print "<a href=\"JavaScript: nextpage($newpage);\"><img SRC=../images/back.gif border=0></a>"; }
else{ print "<img SRC=../images/back1.gif border=0>"; }
my $newpage = $page+1;
if($page < $nPages){ print "<a href=\"JavaScript: nextpage($newpage);\"><img SRC=../images/forward.gif border=0></a>"; }
else{ print "<img SRC=../images/forward1.gif border=0>"; }
print "&nbsp; &nbsp; &nbsp; To page:<INPUT NAME=\"go_to_page\" Size=3 VALUE=1>";
print "<INPUT NAME=button_go TYPE=button VALUE=Go onClick=jump_to_page();>";
print "&nbsp; &nbsp; &nbsp; Find motif:<INPUT NAME=\"search_term\" Size=5>";
print "<INPUT NAME=button_search TYPE=button VALUE=Go onClick=search_motif();>";
print "</FORM>\n";

# Table of motifs
if(!open(INFO, $filenameMotif)){
	print "Motif output file not found!\n"; exit(0);
}
my $count;
my $line="";
my @headers=();
my @headerInd=();
while($line = <INFO>){
	$line =~ s/\n$//;
	my $junk;
	if($line =~ /^Headers/){
		($junk,$junk,@headers) = split("\t",$line);
	}	
	if($line =~ /^>/){
		++$count;
		if($count >= ($page-1)*10+1){ last; }
	}
}
for(my $j=0; $j<@headersGlob; ++$j){
	$headerInd[$j]=-1;
	for(my $i=0; $i<@headers; ++$i){
		if($headersGlob[$j] eq $headers[$i]){ $headerInd[$j]=$i; }
	}
}
print "<TABLE BORDER=0><TR><TD>Name<TD><center>Motif logo\n";
for(my $j=0; $j<@headersGlob; ++$j){
	if($headerInd[$j]>=0){ print "<TD>$headersGlob[$j]"; }
}
#if($headerInd[0]<0){ print "No pattern specified in the motif output"; exit(0); }
my $ID = $outputIDlogo;
for(my $imotif=$start; $imotif<$start+10 && $imotif<=$nMotifs; ++$imotif){
	if(!$line){ last; }
	$line =~ s/\n$//;
	$line =~ s/^\>//;
	my ($motifName,@items) = split(/\t/,$line);
	my $pattern;
	if($headerInd[0]>=0){ $pattern = $items[$headerInd[0]]; }
	else{ $pattern = $items[0]; }
	my $WID = length($pattern)*13;
	if($WID > 600){ $WID = 600; }
	print "<TR><TD>$motifName<TD><a href=\"JavaScript: show_motif('$motifName',$ID);\"><img src=../output/$ID.gif border=0 height=25 width=$WID></a>\n";
	print "<TD><a href=\"JavaScript: show_motif('$motifName',$ID);\">$pattern</a>\n";
	for(my $j=1; $j<@headersGlob; ++$j){
		if($headerInd[$j]<0){ next; }
		print "<TD><center>$items[$headerInd[$j]]\n";
	}
	if($file_motif !~ /^public-/){
		print "<TD><INPUT NAME=\"reverse$imotif\" TYPE=button VALUE=Reverse onClick=\"motif_reverse('$motifName');\">\n";
		print "<TD><INPUT NAME=\"delete$imotif\" TYPE=button VALUE=Delete onClick=\"motif_delete('$motifName');\">\n";
	}
	++$ID;
	while($line = <INFO>){
		if($line =~ /^>/){ last; }
	}
}
close INFO;
print "</TABLE><p>\n";
print "<b>Note:</b> Frequency of clusters of motifs = sum of frequencies of members<p>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  show_motif
#**************************************
{
my $motif = $hashInput{"motif"};
my $seqLogo = $hashInput{"seqLogo"};
my $file_motif = $hashInput{"file_motif"};
my $file_motif_output = $hashInput{"file_motif_output"};

my $filenameMotif = "$PATH_DATA/$loginname-$file_motif";
my $file_motif_source = $file_motif;
if($file_motif =~ /^public-/ && !$hashInput{"data_saved"}){ $filenameMotif = "$PATH_DATA/$file_motif"; }
if($file_motif_output){
	$filenameMotif = "$PATH_OUTPUT/$file_motif_output";
	$file_motif_source = $file_motif_output;
}

if(!$seqLogo){
	# Generate motif-logo
	$seqLogo = get_outputID(1);
	`./motiflogo $filenameMotif $PATH_OUTPUT/$seqLogo.gif -name $motif`;
}

# Print page header
print "<HTML><HEAD><TITLE>Motif $motif</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
my $save_motif = 0;
my @motif_list=();
$save_motif = 1;
if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){
	error_message("Configuration file not found!");
}
while(my $line = <INFO>){
	if(length($line) < 3) { next; }
	my %hash=();
	read_config_line($line,\%hash);
	if($hash{"type_motif"}){ push(@motif_list,$hash{"type_motif"}); }
}
close INFO;
@motif_list = sort @motif_list;
my $list1 = join("\",\"",@motif_list);
if($list1){ $list1 = "\"".$list1."\""; }
print "motif_list = new Array($list1);\n";
print "function motif_check(){\n";
print "	if(!document.cisfinder.pfm_name.value){\n";
print "		alert(\"Enter name for new motif\"); return false;\n";
print "	}\n";
print "	if(document.cisfinder.pfm_end.value-document.cisfinder.pfm_start.value < 6){\n";
print "		alert(\"Selected motif length is too short (<7)\"); return false;\n";
print "	}\n";
print "	if(!document.cisfinder.file_motif_save.selectedIndex){\n";
print "		var filename = prompt(\"Enter file name for motifs\",\"\");\n";
print "		if(!filename){ alert(\"Operation cancelled\"); return false; }\n";
print "		var i;\n";
print "		for(i=0; i<motif_list.length; ++i){\n";
print "			if(filename == motif_list[i]){\n";
print "				alert(\"File name already exist!\"); return false;\n";
print "			}\n";
print "		}\n";
print "		var y=document.createElement('option');\n";
print "		y.text=filename;\n";
print "		y.value=filename;\n";
print "  	try{ document.cisfinder.file_motif_save.add(y,null); }\n";
print "  	catch(ex){ document.cisfinder.file_motif_save.add(y); }\n";
print "  	document.cisfinder.file_motif_save.selectedIndex = document.cisfinder.file_motif_save.length-1;\n";
print "	}\n";
print "	return true;\n";
print "}\n";
print "function motif_save(){\n";
print "	if(!motif_check()){ return false; }\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<h3>Motif $motif from file $file_motif_source</h3>\n";

if(!open(INFO, $filenameMotif)){
	print "Motif file not found!\n"; exit(0);
}
my $line="";
my @headers;
my @headerInd;
my $junk;
while($line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^Headers/){
		($junk,$junk,@headers) = split("\t",$line);
	}	
	if($line =~ /^>$motif\t/){ last; }
}
for(my $j=0; $j<@headersGlob; ++$j){
	$headerInd[$j]=-1;
	for(my $i=0; $i<@headers; ++$i){
		if($headersGlob[$j] eq $headers[$i]){ $headerInd[$j]=$i; }
	}
}
$line =~ s/^\>//;
my ($motifName,@items) = split(/\t/,$line);
my $pattern;
if($headerInd[0]>=0){ $pattern = $items[$headerInd[0]]; }
else{ $pattern = $items[0]; }
print "<TABLE BORDER=1><TR>\n";
for(my $j=1; $j<@headersGlob; ++$j){
	if($headerInd[$j]>=0){ print "<TD width=80><center>$headersGlob[$j]"; }
}
print "<TR>\n";
if(@headers){
	for(my $j=1; $j<@headersGlob; ++$j){
		if($headerInd[$j]<0){ next; }
		print "<TD><center>$items[$headerInd[$j]]\n";
	}
}else{
	for(my $j=1; $j<@items; ++$j){
		print "<TD><center>$items[$j]\n";
	}
}
print "</TABLE>\n";

print "<img src=../output/$seqLogo.gif border=0></a>\n";
print "<TABLE BORDER=0><TR>\n";
my(@ch) = split(//,$pattern);
my $len = @ch;
for(my $i=0; $i<@ch; ++$i){ print "<TD width=23><center>$ch[$i]\n"; }
print "<TR>\n";
for(my $i=1; $i<=@ch; ++$i){ print "<TD><center>$i\n"; }
print "</TABLE>\n";

print "<h3>Position Frequency Matrix (PFM)</h3>\n";
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<INPUT NAME=\"motif\" TYPE=\"hidden\" VALUE=\"$motif\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"save_PFM\">\n";
if($file_motif){
	print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=\"$file_motif\">\n";
}
if($file_motif_output){
	print "<INPUT NAME=\"file_motif_output\" TYPE=\"hidden\" VALUE=\"$file_motif_output\">\n";
}
print "<TABLE BORDER=0>\n";
print "<TR><TD><TD><center>Start position<TD><center>End position<TD><center>Reverse\n";
if($save_motif){ print "<TD><center>Name<TD><center>Add PFM to file\n"; }
print "<TR><TD><INPUT NAME=\"save_pfm_button\" TYPE=\"button\" VALUE=\"Copy motif\" onClick=\"motif_save();\">\n";
print "<TD><center><INPUT NAME=\"pfm_start\" SIZE=2 VALUE=1>\n";
print "<TD><center><INPUT NAME=\"pfm_end\" SIZE=2 VALUE=$len>\n";
print "<TD><center><INPUT NAME=\"reverse_PFM\" TYPE=\"checkbox\">\n";
if($save_motif){
	print "<TD><center><INPUT NAME=\"pfm_name\" SIZE=5>\n";
	print "<TD><center><SELECT name=\"file_motif_save\">\n";
	print "<option value=\"\"> New file\n";
	for(my $i=0; $i<@motif_list; ++$i){ 
		print "<option value=\"$motif_list[$i]\"> $motif_list[$i]";
	}
	print "</select>\n";
}
print "</TABLE>\n";
print "<TABLE BORDER=1><TR><TD><b><font color=blue>Position<TD width=70><center><b><font color=blue>A<TD width=70><center><b><font color=blue>C<TD width=70><center><b><font color=blue>G<TD width=70><center><b><font color=blue>T\n";
while($line = <INFO>){
	$line =~ s/\n$//;
	if(!$line){ last; }
	my @items = split("\t",$line);
	$items[0]++;
	print "<TR>";
	for(my $i=0; $i<5; ++$i){
		print "<TD><center>$items[$i]\n";
	}
}
close INFO;
print "</TABLE>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  search_motif_name
#**************************************
{
my $file_motif = $hashInput{"file_motif1"};
my $file_motif_output = $hashInput{"file_motif_output"};
my $file_motif_elem = $hashInput{"file_motif_elem"};

if(!$file_motif){ $file_motif = $hashInput{"file_motif"}; }
my $filenameMotif = "$PATH_DATA/$loginname-$file_motif";
if($file_motif =~ /^public-/ || $hashInput{"data_saved"}){ $filenameMotif = "$PATH_DATA/$file_motif"; }
if($file_motif_output){
	$filenameMotif = "$PATH_OUTPUT/$file_motif_output";
}elsif(!$file_motif){
	print "No motif file!"; exit(0);
}
my $search_term = $hashInput{"search_term"};
if(!$search_term){
	error_message("No search_term!");
}
my @motifList=();
if(!open(INFO, $filenameMotif)){
	error_message("Motif file $filenameMotif not found!");
}
my $count=0;
while(my $line = <INFO>){
	if($line =~ /^>/){
		++$count;
		my ($name,$junk) = split(/\t/,$line);
		if($name =~ /$search_term/i){
			push(@motifList,$count);
		}
	}
}
close INFO;
my $nMotifs = @motifList;

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "function cleanAll(){\n";
print "	document.cisfinder.motif_delete.value = \"\";\n";
print "	document.cisfinder.motif_reverse.value = \"\";\n";
print "	document.cisfinder.target = \"\";\n";
print "	document.cisfinder.action.value = \"show_motifs\";\n";
print "}\n";
print "function show_motif(motif,logo){\n";
print "	cleanAll();\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.motif.value = motif;\n";
print "	document.cisfinder.seqLogo.value = logo;\n";
if($file_motif_elem){
	print "	document.cisfinder.action.value = \"show_cluster\";\n";
}else{
	print "	document.cisfinder.action.value = \"show_motif\";\n";
}
print "	document.cisfinder.submit();\n";
print "}\n";
print "function motif_delete(motif){\n";
print "	cleanAll();\n";
print "	document.cisfinder.motif_delete.value = motif;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function motif_reverse(motif){\n";
print "	cleanAll();\n";
print "	document.cisfinder.motif_reverse.value = motif;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();

# Generate motif-logo
my $outputIDlogo = get_outputID($nMotifs);
for(my $i=0; $i<$nMotifs; ++$i){
	my $ID = $outputIDlogo+$i;
	`./motiflogo $filenameMotif $PATH_OUTPUT/$ID.gif -num $motifList[$i]`;
}
if(!$nMotifs){
	print "<H3>Search_term not found!</H3>"; exit(0);
}
print "<h3>Motifs in file <font color=BROWN>$file_motif</font> matring to term '$search_term'</h3>\n";
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"search_motif_name\">\n";
print "<INPUT NAME=\"motif\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=$file_motif>\n";
print "<INPUT NAME=\"motif_delete\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"motif_reverse\" TYPE=\"hidden\" VALUE=\"\">\n";
if($file_motif_output){
	print "<INPUT NAME=\"file_motif_output\" TYPE=\"hidden\" VALUE=\"$file_motif_output\">\n";
}
if($file_motif){
	print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=\"$file_motif\">\n";
}
if($file_motif_elem){
	print "<INPUT NAME=\"file_motif_elem\" TYPE=\"hidden\" VALUE=\"$file_motif_elem\">\n";
}
print "<INPUT NAME=\"seqLogo\" TYPE=\"hidden\" VALUE=\"\">\n";
print "</FORM>\n";

# Table of motifs
if(!open(INFO, $filenameMotif)){
	print "Motif output file not found!\n"; exit(0);
}
my $line="";
my @headers=();
my @headerInd=();
while($line = <INFO>){
	$line =~ s/\n$//;
	my $junk;
	if($line =~ /^Headers/){
		($junk,$junk,@headers) = split("\t",$line);
	}	
	if($line =~ /^>/){
		last;
	}
}
for(my $j=0; $j<@headersGlob; ++$j){
	$headerInd[$j]=-1;
	for(my $i=0; $i<@headers; ++$i){
		if($headersGlob[$j] eq $headers[$i]){ $headerInd[$j]=$i; }
	}
}
print "<TABLE BORDER=0><TR><TD>Name<TD><center>Motif logo\n";
for(my $j=0; $j<@headersGlob; ++$j){
	if($headerInd[$j]>=0){ print "<TD>$headersGlob[$j]"; }
}

my $ID = $outputIDlogo;
$count = 1;
for(my $i=0; $i<$nMotifs; ++$i){
	my $ID = $outputIDlogo+$i;
	if(!$line){ last; }
	my $WID;
	my ($motifName,@items);
	my $nLines;
	while($count <= $motifList[$i]){
		$line =~ s/\n$//;
		$line =~ s/^\>//;
		$WID = 0;
		($motifName,@items) = split(/\t/,$line);
		$nLines = 0;
		while($line = <INFO>){
			if($line =~ /^>/){ last; }
			if(length($line) > 1){ ++$nLines; }
		}
		++$count;
	}
	my $pattern;
	if($headerInd[0]>=0){
		$pattern = $items[$headerInd[0]];
		$WID = length($pattern)*13;
	}
	else{
		$pattern = $items[0];
		$WID = $nLines*13;
	}
	if($WID > 600){ $WID = 600; }
	print "<TR><TD>$motifName<TD><a href=\"JavaScript: show_motif('$motifName',$ID);\"><img src=../output/$ID.gif border=0 height=25 width=$WID></a>\n";
	print "<TD><a href=\"JavaScript: show_motif('$motifName',$ID);\">$pattern</a>\n";
	for(my $j=1; $j<@headersGlob; ++$j){
		if($headerInd[$j]<0){ next; }
		print "<TD><center>$items[$headerInd[$j]]\n";
	}
	if($file_motif !~ /^public-/){
		print "<TD><INPUT NAME=\"reverse$i\" TYPE=button VALUE=Reverse onClick=\"motif_reverse('$motifName');\">\n";
		print "<TD><INPUT NAME=\"delete$i\" TYPE=button VALUE=Delete onClick=\"motif_delete('$motifName');\">\n";
	}
	++$ID;
}
close INFO;
print "</TABLE>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  show_search_results
#**************************************
{
my %hashConfig=%hashInput;
if(!$hashInput{"data_saved"}){
	if(!open(INFO, "$PATH_INFO/$loginname-config.txt")){ error_message("Profile not found"); }
	while(my $line = <INFO>){
		if(length($line) < 3) { next; }
		my %hash=();
		read_config_line($line,\%hash);
		if($hashInput{"mainPage"}){
			if($hash{"mainPage"}){ %hashConfig = %hash; last; }
		}else{
			if($dataname && $dataname eq $hash{"data"}){ %hashConfig = %hash; last; }
		}
	}
	close INFO;
	if($hashInput{"mainPage"}){ $dataname=""; }
}
my $file_search = $hashConfig{"file_search"};
my $interval = $hashConfig{"interval"};
if(!$interval){ $interval=100; }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<h3>Specify sequence interval for frequency distribution of motifs</h3>\n";

print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<h3>File with motif search results</h3>\n";
print "<TABLE border=1>\n";
print "<tr><td><b>File type<td><b>File name\n";
print "<tr><td><b>Search file<td><font color=brown>$file_search\n";
print "</TABLE><p>\n";

print "<b>Parameters</b>\n";
my @interval_list = (5,10,20,50,100,200,500,1000,2000);
print "<TABLE BORDER=0>\n";
print "<TR><TD><select name=\"interval\">\n";
for(my $i=0; $i<@interval_list; ++$i){ 
	print "<option value=$interval_list[$i]";
	if($interval==$interval_list[$i]){ print " selected"; }
	print "> $interval_list[$i]\n";
}
print "</select><TD> Interval (bp) for frequency distribution\n";
print "</TABLE>\n";
print "<TR><TD><INPUT NAME=\"submit\" TYPE=submit VALUE=Continue>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"show_search_results1\">\n";
print "<INPUT NAME=\"file_search\" TYPE=\"hidden\" VALUE=\"$file_search\">\n";
if($hashInput{"data_saved"}){
	print "<INPUT NAME=\"data_saved\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
}
if($hashInput{"mainPage"}){
	print "<INPUT NAME=\"mainPage\" TYPE=\"hidden\" VALUE=\"mainPage\">\n";
}
print "</FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  show_search_results1
#**************************************
{
my $file_search = $hashInput{"file_search"};
my $interval = $hashInput{"interval"};

my $fileSearch = $file_search;
if($file_search !~ /^public-/ && !$hashInput{"data_saved"}){ $fileSearch = "$loginname-$file_search"; }
#Read motif file name and fasta file name from the search results file
my($file_motif,$file_fasta);
open(INFO, "$PATH_DATA/$fileSearch");
for(my $i=0; $i<2; ++$i){
	my $line = <INFO>;
	my %hash=();
	if($line =~ /^Parameters:/){
		read_config_line($line,\%hash);
		$file_fasta=$hash{"file_fasta"};
		$file_motif=$hash{"file_motif"};
		$file_fasta =~ s/^$loginname-//;
		$file_motif =~ s/^$loginname-//;
	}
}
close INFO;

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "function show_motifs(){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"show_search_motifs\";\n";
print "	document.cisfinder.search.value = 0;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function show_sequences(){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"show_search_sequences\";\n";
print "	document.cisfinder.search.value = 0;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function search_sequences(){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"show_search_sequences\";\n";
print "	document.cisfinder.search.value = 1;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function get_info(iTF){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	var motif, ind;\n";
print "	if(iTF==1){ ind=document.cisfinder.TF1.selectedIndex; motif=document.cisfinder.TF1.options[ind].value; }\n";
print "	else if(iTF==2){ ind=document.cisfinder.TF2.selectedIndex; motif=document.cisfinder.TF2.options[ind].value; }\n";
print "	else if(iTF==3){ ind=document.cisfinder.TF3.selectedIndex; motif=document.cisfinder.TF3.options[ind].value; }\n";
print "	if(ind==0){\n";
print "		alert(\"Motif is not selected\");\n";
print "		return(false);\n";
print "	}\n";
print "	document.cisfinder.motif.value = motif;\n";
print "	document.cisfinder.action.value = \"show_motif\";\n";
print "	document.cisfinder.search.value = 1;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<h3>Motif search results: $file_search</h3>\n";

print "Fasta file = <font color=brown>$file_fasta</font><br>\n";
print "Motif file = <font color=brown>$file_motif</font><p>\n";

my $outputIDfreq = get_outputID(3);
my $outputIDabund = $outputIDfreq+1;
my $outputIDsearch = $outputIDfreq+2;
my @command=(
	"./patternDistrib",
	"-i $PATH_DATA/$fileSearch",
	"-f $PATH_OUTPUT/$outputIDfreq.txt",
	"-a $PATH_OUTPUT/$outputIDabund.txt",
	"-int $interval",
);
print "<b>Making frequency and abundance tables...</b><br>\n";
my $string = join(' ',@command);
#print "$string\n";
my $response = `$string`;
system("cp $PATH_DATA/$fileSearch $PATH_OUTPUT/$outputIDsearch.txt");
print "<pre>$response</pre>\n";
if($response =~ /ERROR/i){
	print "<h3>Errors in program run<h3>Output was not generated<p>\n";
	exit(0);
}
my @items=split(/ +/,$response);
my $nSeq = $items[3];
my $nMotifs = $items[5];
my $maxPos = $items[7];
if(!$nMotifs){
	print "<b>No motifs in the file!</b><p>\n";
	exit(0);
}

print "<b>Output in text format:</b><br>\n";
print "<a href=../output/$outputIDsearch.txt target=_BLANK1621>Search results</a><br>\n";
print "<a href=../output/$outputIDfreq.txt target=_BLANK1634>Frequency table</a><br>\n";
print "<a href=../output/$outputIDabund.txt target=_BLANK1623>Abundance table</a><p>\n";

# Display results
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"show_search_motifs\">\n";
print "<INPUT NAME=\"file_search\" TYPE=\"hidden\" VALUE=$file_search>\n";
print "<INPUT NAME=\"file_freq\" TYPE=\"hidden\" VALUE=$outputIDfreq.txt>\n";
print "<INPUT NAME=\"file_abund\" TYPE=\"hidden\" VALUE=$outputIDabund.txt>\n";
print "<INPUT NAME=\"file_fasta\" TYPE=\"hidden\" VALUE=$file_fasta>\n";
print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=$file_motif>\n";
print "<INPUT NAME=\"interval\" TYPE=\"hidden\" VALUE=$interval>\n";
print "<INPUT NAME=\"motif\" TYPE=\"hidden\">\n";
print "<INPUT NAME=\"nMotifs\" TYPE=\"hidden\" VALUE=$nMotifs>\n";
print "<INPUT NAME=\"search\" TYPE=\"hidden\" VALUE=0>\n";
print "<INPUT NAME=\"show_frequency\" TYPE=\"hidden\" VALUE=1>\n";
print "<INPUT NAME=\"show_search_motif_button\" TYPE=button VALUE=\"Show all motifs\" onClick=\"show_motifs();\">\n";
print "<INPUT NAME=\"show_search_sequence_button\" TYPE=button VALUE=\"Show all sequences\" onClick=\"show_sequences();\"><p>\n";

print "<HR NOSHADE COLOR=BLUE></HR>\n";

print "<h3><INPUT NAME=\"search_button\" TYPE=button VALUE=\"Search\" onClick=\"search_sequences();\"> &nbsp; &nbsp; &nbsp; &nbsp; Find sequences</h3>\n";
print "<table border=0>\n";
print "<tr><td><b>Name search:</b><td><input size=20 name=search_term>\n";
print "<select name=field>\n";
print "	<option value=default>----- Select -----\n";
print "	<option value=sequence>Sequence name\n";
print "	<option value=symbol>Gene Symbol\n";
print "</select>\n";
print "<tr><td><b>Distance from TSS:</b>\n";
print "<td><select name=distanceTSS><option value=0>No limit<option value=100>100<option value=200>200<option value=500>500<option value=1000>1000<option value=2000>2000<option value=5000>5000<option value=10000>10000<option value=20000>20000<option value=50000>50000<option value=100000>100000<option value=200000>200000\n";
print "</select> (If annotations uploaded)\n";
print "</table>\n";

# Make motif selection menu
my $fileMotif = $file_motif;
if($file_motif !~ /^public-/ && !$hashInput{"data_saved"}){ $fileMotif = "$loginname-$file_motif"; }
my @motifs=();
if(!open(INFO, "$PATH_DATA/$fileMotif")){
	print "Motif file not found!\n"; exit(0);
}
while(my $line = <INFO>){
	if($line =~ />/){
		$line =~ s/\n$//;
		$line =~ s/^>//;
		my($name,$junk)=split(/\t/,$line);
		push(@motifs,$name);
	}
}
close INFO;
my $select_motifs="<OPTION VALUE=\"\">----- Select motif -----\n";
foreach my $name (sort @motifs){
	$select_motifs .= "<OPTION VALUE=$name>$name\n";
}

print "<table><tr><td>No.<td>Motif<td>From<td>To<td>N hits<td>Conser-<br>vation %<td>Score<td>Strand\n";
for(my $iTF=1; $iTF<=3; ++$iTF){
	print "<tr><td>$iTF.<td><SELECT NAME=TF$iTF>$select_motifs\n";
	print "</SELECT>\n";
	print "<td><input name=startTF$iTF size=7 value=0>\n";
	print "<td><input name=endTF$iTF size=7 value=1000000>\n";
	print "<td>N&ge;<select name=nhitsTF$iTF><option value=1>1<option value=2>2<option value=3>3<option value=5>5<option value=10>10\n";
	print "<td><SELECT name=conservTF$iTF>\n";
	for(my $ic=0; $ic<10; ++$ic){
		my $cons = $ic*10;
		print "<OPTION VALUE=$cons>$cons";
	}
	print "</SELECT>\n";
	print "<td><SELECT name=scoreTF1>\n";
	for(my $is=4; $is<12; ++$is){
		print "<OPTION VALUE=$is>$is";
	}
	print "</SELECT>\n";
	print "<td><SELECT name=strandTF1>\n";
	print "<OPTION VALUE=0>N/A<OPTION VALUE=1>+<OPTION VALUE=-1>-\n";
	print "</SELECT>\n";
	print "<td><input type=button onClick=get_info($iTF); value=\"Motif Info\">\n";
}
print "</TABLE>\n";
print "<b>Combine motifs:</b> <SELECT NAME=combineMotifs>\n";
print "<OPTION VALUE=and>And\n";
print "<OPTION VALUE=or>Or\n";
print "<OPTION VALUE=dist>Distance\n";
print "</SELECT>\n";
print "&nbsp; &nbsp;Max distance: <SELECT NAME=distance><OPTION VALUE=10>10<OPTION VALUE=20>20<OPTION VALUE=50 SELECTED>50<OPTION VALUE=100>100<OPTION VALUE=200>200\n";
print "</SELECT>\n";
print "</FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  show_search_motifs
#**************************************
{
my $file_search = $hashInput{"file_search"};
my $file_freq = $hashInput{"file_freq"};
my $file_fasta = $hashInput{"file_fasta"};
my $file_motif = $hashInput{"file_motif"};
my $show_frequency = $hashInput{"show_frequency"};
my $file_search_output = $hashInput{"file_search_output"};

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "function show_motif(motif){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.motif.value = motif;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<h3>Motifs from file $file_motif in sequences $file_fasta</h3>\n";
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"show_search_motif\">\n";
print "<INPUT NAME=\"motif\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=$file_motif>\n";
print "<INPUT NAME=\"file_fasta\" TYPE=\"hidden\" VALUE=$file_fasta>\n";
print "<INPUT NAME=\"file_search\" TYPE=\"hidden\" VALUE=$file_search>\n";
print "<INPUT NAME=\"file_search_output\" TYPE=\"hidden\" VALUE=$file_search_output>\n";
print "<INPUT NAME=\"file_freq\" TYPE=\"hidden\" VALUE=$file_freq>\n";
print "<INPUT NAME=\"show_frequency\" TYPE=\"hidden\" VALUE=$show_frequency>\n";
print "</FORM>\n";
# Table of motifs
my %hashMotif=();
if(!open(INFO, "$PATH_OUTPUT/$file_freq")){
	print "Frequency file not found!\n"; exit(0);
}
my $line = <INFO>;
while($line = <INFO>){
	$line =~ s/\n$//;
	my($motifName,$strand,@values) = split(/\t/,$line);
	if($strand eq "+" || $strand eq "-"){
		foreach my $x (@values){ $hashMotif{$motifName} += $x; };
	}
}
close INFO;
my @motifNames = sort keys %hashMotif;
print "<TABLE BORDER=0><TR><TD>Motif Name<TD><center>Frequency\n";
for(my $j=0; $j<@motifNames; ++$j){
	print "<TR><TD><a href=\"JavaScript: show_motif('$motifNames[$j]');\">$motifNames[$j]</a>\n";
	print "<TD><center>$hashMotif{$motifNames[$j]}\n";
}
print "</TABLE>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  show_search_motif
#**************************************
{
my $motif = $hashInput{"motif"};
my $file_search = $hashInput{"file_search"};
my $file_motif = $hashInput{"file_motif"};
my $file_fasta = $hashInput{"file_fasta"};
my $file_freq = $hashInput{"file_freq"};
my $show_frequency = $hashInput{"show_frequency"};
my $file_search_output = $hashInput{"file_search_output"};

my $fileSearch = "$PATH_DATA/$file_search";
if($file_search !~ /^public-/ && !$hashInput{"data_saved"}){ $fileSearch = "$PATH_DATA/$loginname-$file_search"; }
if($file_search_output){ $fileSearch = "$PATH_OUTPUT/$file_search_output"; }

my $outputIDall = get_outputID(3);
my $outputIDstrand = $outputIDall+1;
my $outputIDsites = $outputIDall+2;
if($show_frequency){
	plot_frequency_distribution("$PATH_OUTPUT/$file_freq",$motif,$outputIDall);
}
# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "function show_motif(){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();

print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"data\" TYPE=\"hidden\" VALUE=\"$dataname\">\n";
print "<INPUT NAME=\"motif\" TYPE=\"hidden\" VALUE=\"$motif\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"show_motif\">\n";
print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=\"$file_motif\">\n";
print "</FORM>\n";

print "<h3>Motif <a href=\"JavaScript: show_motif();\">$motif</a> (from $file_motif) in sequences from $file_fasta</h3>\n";
print "<a href=../output/$outputIDsites.txt target=_blank11>Full table of matches</a> - in a text format<p>\n";
if($show_frequency){
	print "<b>Frequency distribution in all sequences (100,50,25, and 12% quantiles by score) </b><br>\n";
	print "<img src=../output/$outputIDall.gif><p>\n";
	print "<b>Frequency distribution by strand</b><br>\n";
	print "<img src=../output/$outputIDstrand.gif><p>\n";
}
my @hits=();
if(!open(INFO, $fileSearch)){
	print "Search file not found!\n"; exit(0);
}
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^Headers|^Parameters/){ next; }
	my($motifName,$seqName,$strand,$len,$pos,$score,$match,$cons) = split(/\t/,$line);
	if($motifName eq $motif){
		push(@hits,[$seqName,$strand,$len,$pos,$score,$match,$cons]);
	}
}
close INFO;
my $fileConserv = "$PATH_DATA/$file_fasta";
if($file_fasta !~ /^public-/ && !$hashInput{"data_saved"}){ $fileConserv="$PATH_DATA/$loginname-$file_fasta"; }
$fileConserv =~ s/\.fa$/.cons/;
my $consExist = file_exist($fileConserv);
my @hits = sort {$b->[4]<=>$a->[4]} @hits;
print "<TABLE BORDER=0><TR><TD><b><font color=blue>No.<TD><b><font color=blue>Motif Name<TD><b><font color=blue>Sequence Name<TD><b><font color=blue>Strand<TD><b><font color=blue>Length\n";
print "<TD><b><font color=blue>Start<TD width=80><center><b><font color=blue>Score<TD><b><font color=blue>p-value<TD><b><font color=blue>Pattern\n";
if($consExist){ print "<TD><b><font color=blue>Conservation\n"; }
my $log10 = log(10.0);
open(OUT, ">$PATH_OUTPUT/$outputIDsites.txt");
for(my $i=0; $i<@hits; ++$i){
	my $ref = $hits[$i];
	my ($seqName,$strand,$len,$pos,$score,$match,$cons) = @$ref;
	my $p = int(100000*exp(-$score*$log10))/100000;
	my $count = $i+1;
	if($i<2000){
		print "<TR><TD>$count<TD>$motif<TD><center>$seqName<TD><center>$strand<TD><center>$len<TD><center>$pos<TD><center>$score<TD><center>$p<TD><font face=courier>$match\n";
		if($consExist){ print "<TD><center>$cons\n"; }
	}
	print OUT "$count\t$motif\t$seqName\t$strand\t$len\t$pos\t$score\t$p\t$match";
	if($consExist){ print OUT "\t$cons"; }
	print OUT "\n";
}
close OUT;
print "</TABLE>\n";
if(@hits>2000){
	print "The number of records exceeds 2000<br>\n";
	print "<a href=../output/$outputIDsites.txt target=_blank11>Full table of matches</a> in a text format\n";
}
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  show_search_sequences
#**************************************
{
my $file_search = $hashInput{"file_search"};
my $file_motif = $hashInput{"file_motif"};
my $file_fasta = $hashInput{"file_fasta"};
my $page = $hashInput{"page"};
my $nPages = $hashInput{"nPages"};
my $nSequences = $hashInput{"nSequences"};
my $file_search_output = $hashInput{"file_search_output"};
my $file_abund = $hashInput{"file_abund"};
my $file_sequence_table = $hashInput{"file_sequence_table"};

#Make table of sequences
my %hashSeq=();
if(!open(INFO, "$PATH_OUTPUT/$file_abund")){
	error_message("Abundance file not found!");
}
my $line = <INFO>;
while($line = <INFO>){
	$line =~ s/\n$//;
	my($name,@items) = split(/\t/,$line);
	my $N = 0;
	foreach my $k (@items){ $N += $k; }
	$hashSeq{$name}->[0] = $N;
}
close INFO;

my @seqNames = sort(keys %hashSeq);
if(!$page || !$nSequences || !$file_sequence_table){
	my $fileFasta = $file_fasta;
	if($file_fasta !~ /^public-/ && !$hashInput{"data_saved"}){
		$fileFasta = "$loginname-$file_fasta";
	}
	my $fileCoordinates = $fileFasta;
	if($fileFasta =~ /\.fa$/){ $fileCoordinates =~ s/\.fa$/.coord/; }
	else{ $fileCoordinates .= ".coord"; }
	my $fileAttributes = $fileFasta;
	if($fileFasta =~ /\.fa$/){ $fileAttributes =~ s/\.fa$/.attr/; }
	else{ $fileAttributes .= ".attr"; }
	my $header_line = "No\tSequence Name\tN motifs";
	if(open(INFO, "$PATH_DATA/$fileCoordinates")){
		$header_line .= "\tChr\tStrand\tLength\tStart";
		my $found=0;
		while(my $line = <INFO>){
			$line =~ s/\n$//;
			if($line =~ /^Headers|^Parameters|genome/){ next; }
			my ($name,$chr,$strand,$len,$pos) = split(/\t/,$line);
			if(defined($hashSeq{$name})){
				push(@{$hashSeq{$name}},$chr,$strand,$len,$pos);
				$found=1;
			}
		}
		close INFO;
		if(!$found){ print "<b>Note: genome location not found</b><p>\n"; }
	}
	if(open(INFO, "$PATH_DATA/$fileAttributes")){
		my $found=0;
		my $line = <INFO>;
		$line =~ s/\n$//;
		my($keyword,$junk,@headers) = split(/\t/,$line);
		if($keyword !~ /^Headers/){ error_message("Headers missing in attribute file"); }
		$header_line .= "\t".join("\t",@headers);
		while(my $line = <INFO>){
			$line =~ s/\n$//;
			my($name,@values) = split(/\t/,$line);
			if(defined($hashSeq{$name})){
				push(@{$hashSeq{$name}},@values);
				$found=1;
			}
		}
		close INFO;
		if(!$found){ print "<b>Note: atributes of sequences not found</b><p>\n"; }
	}
	if($hashInput{"search"}){
		advanced_search(\@seqNames,\%hashSeq,$header_line);
	}
	my $outputIDtable = get_outputID(1);
	$file_sequence_table = "$outputIDtable.txt";
	open(OUT, ">$PATH_OUTPUT/$file_sequence_table");
	print OUT $header_line."\n";
	for(my $j=0; $j<@seqNames; ++$j){
		my $seqName = $seqNames[$j];
		my $no = $j+1;
		print OUT "$no\t$seqName\t".join("\t",@{$hashSeq{$seqName}})."\n";
	}
	close OUT;
	$page=1;
	$nSequences= @seqNames;
	$nPages = int(($nSequences+99)/100);
}

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "var Npages = $nPages;\n";
print "function nextpage(page){\n";
print "	document.cisfinder.page.value = page;\n";
print "	document.cisfinder.action.value = \"show_search_sequences\";\n";
print "	document.cisfinder.target = \"\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function jump_to_page(){\n";
print "	var page = document.cisfinder.go_to_page.value;\n";
print "	if(isNaN(page)){ alert(\"Wrong page number!\"); return false; }\n";
print "	page = Math.ceil(page);\n";
print "	if(page<1 || page>Npages){ alert(\"Wrong page number!\"); return false; }\n";
print "	nextpage(page);\n";
print "}\n";
print "function show_sequence(sequence){\n";
print "	var x = Math.round(Math.random()*1000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"show_search_sequence\";\n";
print "	document.cisfinder.sequence.value = sequence;\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";
print_header();
print "<h3>Sequence list used for motif search (N=$nSequences; Page $page out of $nPages)</h3>\n";
print "Fasta file = <font color=brown>$file_fasta</font><br>\n";
print "Motif file = <font color=brown>$file_motif</font>\n";

print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
my $newpage = $page-1;
if($page > 1){ print "<a href=\"JavaScript: nextpage($newpage);\"><img SRC=../images/back.gif border=0></a>"; }
else{ print "<img SRC=../images/back1.gif border=0>"; }
my $newpage = $page+1;
if($page < $nPages){ print "<a href=\"JavaScript: nextpage($newpage);\"><img SRC=../images/forward.gif border=0></a>"; }
else{ print "<img SRC=../images/forward1.gif border=0>"; }
print "&nbsp; &nbsp; &nbsp; To page:<INPUT NAME=\"go_to_page\" Size=3 VALUE=1>";
print "<INPUT NAME=button_go TYPE=button VALUE=Go onClick=jump_to_page();>";
print "&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; <b>Table as text:</b> <a href=../output/$file_sequence_table target=_BLANK397>table</a><p>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"show_search_sequence\">\n";
print "<INPUT NAME=\"sequence\" TYPE=\"hidden\" VALUE=\"\">\n";
print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=$file_motif>\n";
print "<INPUT NAME=\"file_fasta\" TYPE=\"hidden\" VALUE=$file_fasta>\n";
print "<INPUT NAME=\"file_search\" TYPE=\"hidden\" VALUE=$file_search>\n";
print "<INPUT NAME=\"page\" TYPE=\"hidden\" VALUE=$page>\n";
print "<INPUT NAME=\"nPages\" TYPE=\"hidden\" VALUE=$nPages>\n";
print "<INPUT NAME=\"nSequences\" TYPE=\"hidden\" VALUE=\"$nSequences\">\n";
print "<INPUT NAME=\"file_search_output\" TYPE=\"hidden\" VALUE=\"$file_search_output\">\n";
print "<INPUT NAME=\"file_sequence_table\" TYPE=\"hidden\" VALUE=\"$file_sequence_table\">\n";
my $TF1 = $hashInput{"TF1"};
if($TF1){
	print "<INPUT NAME=\"TF1\" TYPE=\"hidden\" VALUE=\"$TF1\">\n";
	my $TF2 = $hashInput{"TF2"};
	if($TF2){
		print "<INPUT NAME=\"TF2\" TYPE=\"hidden\" VALUE=\"$TF2\">\n";
		my $TF3 = $hashInput{"TF3"};
		if($TF3){ print "<INPUT NAME=\"TF3\" TYPE=\"hidden\" VALUE=\"$TF3\">\n"; }
	}
}
print "</FORM>\n";

if(!open(INFO, "$PATH_OUTPUT/$file_sequence_table")){
	print "Table file not found!\n"; exit(0);
}
my $line = <INFO>;
$line =~ s/\n$//;
my @items = split(/\t/,$line);
print "<TABLE><TR><TD><B>".join("<TD><center><b>",@items)."\n";
my $count = 0;
while(my $line = <INFO>){
	my $page1 = int($count/100);
	++$count;
	if($page1 < $page-1){ next; }
	if($page1 > $page-1){ last; }
	$line =~ s/\n$//;
	my ($no,$name,@items) = split(/\t/,$line);
	print "<TR><TD>$no<TD><a href=\"JavaScript: show_sequence('$name');\">$name</a><TD><center>";
	print join("<td><center>",@items)."\n";
}
print "</TABLE>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  show_search_sequence
#**************************************
{
my $dataname = $hashInput{"data"};
my $sequence = $hashInput{"sequence"};
my $file_search = $hashInput{"file_search"};
my $file_motif = $hashInput{"file_motif"};
my $file_fasta = $hashInput{"file_fasta"};
my $file_search_output = $hashInput{"file_search_output"};

my $fileSearch = "$PATH_DATA/$file_search";
if($file_search !~ /^public-/ && !$hashInput{"data_saved"}){ $fileSearch = "$PATH_DATA/$loginname-$file_search"; }
if($file_search_output){ $fileSearch = "$PATH_OUTPUT/$file_search_output"; }

my @hits=();
my %hashMotif=();
if(!open(INFO, $fileSearch)){
	error_message("Search file not found!");
}
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^Headers|^Parameters/){ next; }
	if($line =~ /^NONE\t/){ next; }
	my($motifName,$seqName,$strand,$len,$pos,$score,$match,$cons) = split(/\t/,$line);
	if($seqName eq $sequence){
		push(@hits,[$motifName,$strand,$len,$pos,$score,$match,$cons]);
		$hashMotif{$motifName}=1;
	}
}
close INFO;
my @hits = sort {$a->[3]<=>$b->[3]} @hits;
my @motif_list = sort keys %hashMotif;

my $fileFasta = $file_fasta;
if($file_fasta !~ /^public-/ && !$hashInput{"data_saved"}){
	$fileFasta = "$loginname-$file_fasta";
}
my $fileCoordinates = $fileFasta;
$fileCoordinates =~ s/\.fa$/.coord/;
my $fileAttributes = $fileFasta;
$fileAttributes =~ s/\.fa$/.attr/;
my $fileConservation = $fileFasta;
$fileConservation =~ s/\.fa$/.cons/;


# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\">\n";
print "<!--\n";
print "function show_motif(motif){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.motif.value = motif;\n";
print "	document.cisfinder.action.value = \"show_motif\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "function show_sequence(){\n";
print "	var x = Math.round(Math.random()*10000);\n";
print "	document.cisfinder.target = \"_BLANK\"+x;\n";
print "	document.cisfinder.action.value = \"show_search_sequence_seq\";\n";
print "	document.cisfinder.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT>\n";

print_header();
print "<h3>Motifs found in sequence $sequence</h3>\n";
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"show_search_sequence_seq\">\n";
print "<INPUT NAME=\"sequence\" TYPE=\"hidden\" VALUE=\"$sequence\">\n";
print "<INPUT NAME=\"motif\" TYPE=\"hidden\">\n";
print "<INPUT NAME=\"file_motif\" TYPE=\"hidden\" VALUE=$file_motif>\n";
print "<INPUT NAME=\"file_fasta\" TYPE=\"hidden\" VALUE=$file_fasta>\n";
print "<INPUT NAME=\"file_search\" TYPE=\"hidden\" VALUE=$file_search>\n";
print "<INPUT NAME=\"file_search_output\" TYPE=\"hidden\" VALUE=$file_search_output>\n";
print "<TABLE BORDER=0>\n";
for(my $iTF=1; $iTF<=3; ++$iTF){
	print "<TR><TD>Motif #$iTF\n";
	my $TF = $hashInput{"TF$iTF"};
	print "<TD><SELECT name=\"motif$iTF\">\n";
	print "<option value=\"\"> ------------------------------- select -------------------------------\n";
	foreach my $motif (@motif_list){ 
		print "<option value=\"$motif\"";
		if($motif eq $TF){ print " SELECTED"; }
		print ">$motif\n";
	}
	print "</select>\n";
}
print "<TR><TD>Show sequence:<TD><INPUT NAME=\"show_sequence_button\" TYPE=button VALUE=\"Show Sequence\" onClick=show_sequence();>\n";
print "</TABLE>\n";
print "</FORM>\n";

my $found = 0;
if(open(INFO, "$PATH_DATA/$fileCoordinates")){
	my ($genome,$chr,$strand,$len,$pos,$junk);
	while(my $line = <INFO>){
		$line =~ s/\n$//;
		if($line =~ /^Headers|^Parameters|genome/){
			if($line =~ /genome/i){
				$line =~ s/.*genome\s*=\s*//i;
				($genome,$junk) = split(/\t/,$line);
			}
		}elsif($line =~ /^$sequence\t/){
			($junk,$chr,$strand,$len,$pos) = split(/\t/,$line);
			$found = 1;
			last;
		}
	}
	close INFO;
	if(!$genome){ print "<b>Note: genome is not specified in the file with sequence coordinates</b><p>\n"; }
	elsif(!$found){ print "<b>Note: coordinates of sequence $sequence not found</b><p>\n"; }
	elsif($chr !~ /^chr/ || $len<5 || $strand !~ /\+|-/){ print "<b>Note: error in coordinates</b><p>\n"; }
	else{
		my $end = $pos+$len;
		my $center = int($pos+$len/2);
		print "<b>Sequence location</b><br>\n";
		print "<table border=1><tr><td>Genome<td><center>Chr<td>Strand<td><center>Start<td><center>End<TD><center><center>Genome<br>Browser<TD><center>Browser of<br>regulatory regions\n"; 
		print "<tr><td><center>$genome<td><center>$chr<td><center>$strand<td><center>$pos<td><center>$end\n"; 
		print "<td><center><a href=\"http://genome.ucsc.edu/cgi-bin/hgTracks?position=$chr:$pos-$end&db=$genome\" target=_BLANK45>UCSC</a><TD>\n";
		if($genome eq "mm9"){
			print "<center><a href=/geneindex/$genome/bin/giReg.cgi?userSeqname=$sequence&chr=$chr&position=$center&userStrand=%2B&Qsize=$len&userBlockSizes=$len,&userTstarts=$pos target=_BLANK64>CisView</a>\n";
		}
		print "</table><p>\n";
	}
}
if(open(INFO, "$PATH_DATA/$fileAttributes")){
	my @values=();
	my $line = <INFO>;
	$line =~ s/\n$//;
	my($junk,$junk1,@headers) = split(/\t/,$line);
	while(my $line = <INFO>){
		if($line =~ /^$sequence\t/){
			$line =~ s/\n$//;
			my($junk,@items) = split(/\t/,$line);
			@values = @items;
			$found = 1;
			last;
		}
	}
	close INFO;
	if(!$found){ print "<b>Note: atributes of sequence $sequence not found</b><p>\n"; }
	else{
		print "<b>Sequence attributes</b><br>\n";
		print "<table border=1><tr><td><center>".join("<td><center>",@headers)."\n";
		print "<tr><td><center>".join("<td><center>",@values)."\n";
		print "</table><p>\n";
	}
}
my $consExist = file_exist("$PATH_DATA/$fileConservation");

# Table of motifs
print "<h4>Motif scanning results</h4>\n";
print "<TABLE BORDER=0><TR><TD><b><font color=blue>No.<TD><b><font color=blue>Motif Name<TD><b><font color=blue>Sequence Name<TD><b><font color=blue>Strand<TD><b><font color=blue>Length\n";
print "<TD><b><font color=blue>Start<TD width=80><center><b><font color=blue>Score<TD><b><font color=blue>p-value<TD><b><font color=blue>Pattern\n";
if($consExist){ print "<TD><b><font color=blue>Conservation\n"; }
my $count=1;
my $icons = 0;
my $log10 = log(10.0);
foreach my $ref (@hits){
	my ($motifName,$strand,$len,$pos,$score,$match,$cons) = @$ref;
	my $p = int(100000*exp(-$score*$log10))/100000;
	print "<TR><TD>$count";
	print "<TD><a href=\"JavaScript: show_motif('$motifName');\">$motifName</a>";
	print "<TD><center>$sequence<TD><center>$strand<TD><center>$len<TD><center>$pos<TD><center>$score<TD><center>$p<TD><font face=courier>$match\n";
	if($consExist){ print "<TD><center>$cons\n"; }
	++$count;
}
print "</TABLE>\n";
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  show_search_sequence_seq
#**************************************
{
my $sequence = $hashInput{"sequence"};
my $file_search = $hashInput{"file_search"};
my $file_motif = $hashInput{"file_motif"};
my $file_fasta = $hashInput{"file_fasta"};
my $motif1 = $hashInput{"motif1"};
my $motif2 = $hashInput{"motif2"};
my $motif3 = $hashInput{"motif3"};
my $file_search_output = $hashInput{"file_search_output"};

my $fileSearch = "$PATH_DATA/$file_search";
if($file_search !~ /^public-/ && !$hashInput{"data_saved"}){ $fileSearch = "$PATH_DATA/$loginname-$file_search"; }
if($file_search_output){ $fileSearch = "$PATH_OUTPUT/$file_search_output"; }
my $fileFasta = $file_fasta;
if($file_fasta !~ /^public-/ && !$hashInput{"data_saved"}){ $fileFasta = "$loginname-$file_fasta"; }

# Print page header
print "<HTML><HEAD><TITLE>CisFinder results $dataname</TITLE>\n";
print_header();
print "<h3>Sequence $sequence with selected motifs</h3>\n";

if(!open(INFO, "$PATH_DATA/$fileFasta")){
	print "Fasta file $file_fasta not found!<br>\n"; exit(0);
}
while(my $line = <INFO>){
	if($line eq ">$sequence\n"){ last; }
}
my $seq = "";
while(my $line = <INFO>){
	if($line =~ /^>/){ last; }
	$line =~ s/\n$//;
	$seq .= $line;
}
close INFO;
my $length = length($seq);
if(!open(INFO, $fileSearch)){
	print "Search file not found!\n"; exit(0);
}
my @hits=();
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^Headers|^Parameters/){ next; }
	if($line =~ /^NONE\t/){ next; }
	my($motifName,$seqName,$strand,$len,$pos,$score,$match) = split(/\t/,$line);
	if($seqName eq $sequence){
		my $ind = 0;
		if($motifName eq $motif1){ $ind=1; }
		elsif($motifName eq $motif2){ $ind=2; }
		elsif($motifName eq $motif3){ $ind=3; }
		if($ind){	
			push(@hits,[$pos,$len,$ind,$strand,$score]);
		}
	}
}
push(@hits,[$length+1,0,0,1,0]);
close INFO;
my @hits = sort {$a->[0]<=>$b->[0]} @hits;
if($motif1 || $motif2 || $motif3){
	print "<TABLE BORDER=0 CELLPADDING=5><TR><TD><b>Motifs:\n";
	if($motif1){ print "<TD><FONT COLOR=MAGENTA>$motif1"; }
	if($motif2){ print "<TD><FONT COLOR=BLUE>$motif2"; }
	if($motif3){ print "<TD><FONT COLOR=ORANGE>$motif3"; }
	print "</TABLE>\n";
}
my $start = 0;
my $im = 0;
my $count = 0;
my $over=0;
my $Nhits = @hits;
my $length = length($seq);
my $order = int(log($length)/log(10))+1;
print "<font face=courier><font color=green>Pos";
for(my $i=0; $i<$order-1; ++$i){ print "&nbsp;"; }
print "<u>";
for(my $i=0; $i<100; ++$i){ my $k=$i%10; print $k; }
print "</u></font><br>\n";
while(1){
	my $ind = sprintf("%0$order"."d",$count);
	print $ind."&nbsp;&nbsp;";
	my $end = $start+100;
	if($end > $length){ $end=$length; }
	my $x = $hits[$im]->[0];
	$over = 0;
	while($x < $end){
		#print "B $x $end $start $im<br>\n";
		if($x>$start){ print substr($seq,$start,$x-$start); }
		if($hits[$im]->[2]==1){ print "<FONT COLOR=MAGENTA>"; }
		elsif($hits[$im]->[2]==2){ print "<FONT COLOR=BLUE>"; }
		else{ print "<FONT COLOR=ORANGE>"; }
		my $endField = $x + $hits[$im]->[1];
		if($endField > $end){
			$endField = $end;
			$over = 1;
		}
		if($im < $Nhits-1 && $hits[$im+1]->[0] < $endField){
			$endField = $hits[$im+1]->[0];
		}
		if($x < $start){ $x=$start; }
		if($endField > $x){ print substr($seq,$x,$endField-$x); }
		print "</FONT>";
		$start = $endField;
		++$im;
		$x = $hits[$im]->[0];
	}
	if($end>$start){ print substr($seq,$start,$end-$start); }
	print "<br>\n";
	if($end==$length){ last; }
	if($over){ --$im; }
	$start = $end;
	$count += 100;;
}
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#***********************************
sub file_exist
#***********************************
{
	if(open (INFO_TEMP, $_[0])){ close INFO_TEMP; 1; }
	else { 0; }
}

#***********************
sub validate
#***********************
{
my $loginname = shift;
my $passwd = shift;

my $line;
my @abcd;
if($SECURITY<1){ return(0); }
if ($loginname =~ /^guest/){ return(0); }
if (!$loginname){
	error_message("You need to put your name in the form!");
}
if(!open(INFO, "$PATH_INFO/login.txt")){
	error_message("Validation failed!");
}
my $metaTerm = quotemeta($loginname);
while($line = <INFO>){
	if($line =~ /^$metaTerm\t/){ last; };
}
close INFO;
@abcd = split (/\t/, $line);
my $encryptedPsw = $passwd;
for(my $i=0; $i<47; $i++){
	$encryptedPsw = crypt $encryptedPsw, $abcd[1];
}
if($loginname ne $abcd[0] || $encryptedPsw ne $abcd[1]){
	my $text = "Incorrect login name or password\n";
	if($line =~ /^$metaTerm\t/){
		$text .= "<FORM METHOD=POST ACTION=$CGI_ADDRESS/cisfinder.cgi>\n";
		$text .= "To reset your password, enter your email: &nbsp; ";
		$text .= "<INPUT NAME=email_reset SIZE=20>";
		$text .= "<INPUT TYPE=submit NAME=reset_password VALUE=\" Reset password \">\n";
		$text .= "<INPUT TYPE=hidden NAME=loginname VALUE=$loginname></FORM>\n";
	}
	error_message($text,"register");
}
return(0);
}

#***********************
sub reset_password
#***********************
{
my $loginname = shift;

my $register="register";
my $line;
my @abcd;
my $email_reset = $hashInput{"email_reset"};
if(!open(INFO, "$PATH_INFO/login.txt")){
	error_message("Operation faled",$register);
}
if(!open(OUT, ">$PATH_INFO/login-temp.txt")){
	close INFO;
	error_message("Operation faled",$register);
}
my @letters = ('A' .. 'Z', 'a' .. 'z', '0' .. '9', '_');
my $salt = $letters[rand@letters] . $letters[rand@letters];
my $content="";
my $email_address;
while($line = <INFO>){
	chop $line;
	@abcd = split (/\t/, $line);
	if($abcd[0] eq $loginname){
		if($abcd[4] !~ /^\S+@\S+$/ || $abcd[4] ne $email_reset){
			error_message("Wrong email address",$register);
		}
		$email_address = $abcd[4];
		my $new_passwd;
		for(my $i=0; $i<8; $i++){
			$new_passwd .= $letters[rand@letters];
		}
		my $encryptedPsw = $new_passwd;
		for(my $i=0; $i<47; $i++){
			$encryptedPsw = crypt $encryptedPsw, $salt;
		}
		$abcd[1]=$encryptedPsw;
		$content .= "Your password was reset:\nloginname=$loginname passwrd=$new_passwd\n\nAutomated mail. Do not reply!\n";
		$line = join("\t",@abcd);
	}
	print OUT $line."\n";
}
close INFO;
close OUT;
my $text;
if($content){
	`echo "$content" | mailx -s "CisFinder reset" $email_address`;
	$text = "Password was reset and sent to you by email!";
	`cp $PATH_INFO/login-temp.txt $PATH_INFO/login.txt`;
	`rm $PATH_INFO/login-temp.txt -f`;
}else{
	$text .= "Operation failed!\n";
}
terminal_window("<H3>$text</H3>",$register);
return(0);
}

#***********************************
sub read_configuration
#***********************************
{
my $filename = shift;
my $hashIni_ref = shift;

if(!open (INFO_TEMP, $filename)){
	print "content-type: text/html","\n\n";
	error_message("Configuration file not found.");
}
while(my $line=<INFO_TEMP>){
	$line =~ s/\n$//;
	my ($keyword,$value) = split(/=/,$line);
	$hashIni_ref->{$keyword} = $value;
}
close INFO_TEMP;
return;
}

#**************************************
sub substitute_chars
#**************************************
{
my $word = shift;
$word =~ s/%2F/\//g;
$word =~ s/\+/ /g;
$word =~ s/%3A/:/g;
$word =~ s/%7E/\~/g;
$word =~ s/%60/\`/g;
$word =~ s/%21/\!/g;
$word =~ s/%23/\#/g;
$word =~ s/%24/\$/g;
$word =~ s/%5E/\^/g;
$word =~ s/%26/\&/g;
$word =~ s/%2B/\+/g;
$word =~ s/%3D/=/g;
$word =~ s/%22/"/g;
$word =~ s/%27/'/g;
$word =~ s/%3B/;/g;
$word =~ s/%3F/\?/g;
$word =~ s/%5C/\\/g;
$word =~ s/%7C/|/g;
$word =~ s/%3C/</g;
$word =~ s/%3E/>/g;
$word =~ s/%2D/-/g;
$word =~ s/%2C/,/g;
$word =~ s/%40/@/g;
$word =~ s/%25/\%/g;
$word =~ s/%28/\(/g;
$word =~ s/%29/\)/g;
return $word;
}

#**************************************
sub print_registration_form
#**************************************
{
print "<title>Registration Form</title>\n";
print_header();
print "<H2>Registration Form</H2>\n";
print "Please fill ALL the fields.<p>\n";
print "<FORM  METHOD=POST ACTION=$CGI_ADDRESS/login.cgi>\n";
print "First name: <input size=20 name=\"firstname\"><p>\n";
print "Last name: &nbsp;<input size=20 name=\"lastname\"><p>\n";
print "Login name: <input size=20 name=\"loginname\"><p>\n";
print "Password: &nbsp; <input type=\"password\" size=20 name=\"passwd\"> <p>\n";
print "Your e-mail: <INPUT SIZE=40 NAME=\"email\"><p>\n";
print "<INPUT TYPE = SUBMIT NAME=\"register\" VALUE = \"Register\"></FORM><p>\n";
exit(0);
}

#**************************************
sub register_new_user
#**************************************
{
my $firstname = uc $hashInput{"firstname"};
my $lastname = uc $hashInput{"lastname"};
my $new_name = $hashInput{"loginname"};
my $new_passwd = $hashInput{"passwd"};
if($new_name =~ /^guest|^public/i){
	error_message("This is not a valid login name. Select another name or login as \"guest\" without password.","register");
}
my $email = $hashInput{"email"};
if($SECURITY>1 && $new_name ne "administrator"){
	if(!open(INFO, "$PATH_INFO/login.txt")){
		error_message("Administrator with loginname \"administrator\" should register first!","register");
	}
	my $line;
	while($line = <INFO>){
		if($line =~ /^administrator\t/){ last; };
	}
	close INFO;
	if ($line !~ /^administrator\t/) {
		error_message("Administrator should register first!","register");
	}
	error_message("Cannot register without administrator!","register");
}
if (!$firstname || !$lastname || !$new_name || !$new_passwd || !$email){
	error_message("Form not fully filled. Registration failed!","register");
}
if ($email !~ /\S\@\S+\.\S/){
	error_message("Registration failed: E-mail address is incorrect!","register");
}
if(open(INFO, "$PATH_INFO/login.txt")){
	my $line;
	my $metaTerm = quotemeta($new_name);
	while($line = <INFO>){
		if($line =~ /^$metaTerm\t/) { last; };
	}
	close INFO;
	if ($line =~ /^$metaTerm\t/) {
		error_message("Login name $new_name is already in use! Select another.","register");
	}
}
my $date;
$date = `date \'+%y%m%d\'`;
chop($date);
if(!open(INFO, ">>$PATH_INFO/login.txt")){
	error_message("Cannot append  file login.txt","register");
}
my @letters = ('A' .. 'Z', 'a' .. 'z', '0' .. '9', '_', '.');
my $salt = $letters[rand@letters] . $letters[rand@letters];
my $encryptedPsw = $new_passwd;
for(my $i=0; $i<47; $i++){
	$encryptedPsw = crypt $encryptedPsw, $salt;
}
print INFO "$new_name\t$encryptedPsw\t$lastname\t$firstname\t$email\t$date\n";
close (INFO);
$sessionID = start_session($new_name);
$loginname = $new_name;
$passwd = $new_passwd;
`cp $PATH_DATA/default-config.txt $PATH_INFO/$loginname-config.txt`;
terminal_window("<H3>Your registration is successful</H3>","continue");
}

#**************************************
sub  read_config_line
#**************************************
{
my $line = shift;
my $hash_ref = shift;

%$hash_ref=();
$line =~ s/\n$//;
my @items = split(/\t/,$line);
foreach my $item (@items){
	my($key,$value) = split(/=/,$item);
	$hash_ref->{$key}=$value;
}
return;
}

#***********************************
sub plot_line
#***********************************
# variable: 0=x1 1=x2 2=color 3=width 4=y_shift 5=scale 6=direction 7=number1 8=number2
{
	my $y1;
	my $x1 = int($_[0]/$_[5]);
	my $x2 = int($_[1]/$_[5]);
	if($_[3] > 1) { $x2 -= 5*$_[6]; }
	my $ystart = int($_[3]/2);
	for(my $i=0; $i<$_[3]; ++$i){
		my $y2 = $_[4]+$ystart-$i;
		print INFO6 "$LINE $_[2] 1 0 $x1 $y2 $x2 $y2\n";
	}
	if($printNumbers){
		$y1 = $_[4] + 7;
		if($x2 > $x1 && $x2 < $x1+18) { $x2 = $x1+18; }
		elsif($x1 >= $x2 && $x1 < $x2+18) { $x1 = $x2+18; }
		if($_[7] != -1){
			print INFO6 "$TEXT $black $tinyFont $x1 $y1 $_[7]\n";
		}
		if($_[8] != -1){
			print INFO6 "$TEXT $black $tinyFont $x2 $y1 $_[8]\n";
		}
	}
}

#***********************************
sub plot_line1
#***********************************
# variable: 0=x1 1=x2 2=color 3=width 4=y_shift 5=scale 6=direction 7=number1 8=number2 9=script
{
	my $y1;
	my $x1 = int($_[0]/$_[5]);
	my $x2 = int($_[1]/$_[5]);
	if($_[3] > 1) { $x2 -= 5*$_[6]; }
	my $ystart = int($_[3]/2);
	my $ref = $_[9];
	for(my $i=0; $i<$_[3]; ++$i){
		my $y2 = $_[4]+$ystart-$i;
		$$ref .= "$LINE $_[2] 1 0 $x1 $y2 $x2 $y2\n";
	}
	if($printNumbers){
		$y1 = $_[4] + 7;
		if($x2 > $x1 && $x2 < $x1+18) { $x2 = $x1+18; }
		elsif($x1 >= $x2 && $x1 < $x2+18) { $x1 = $x2+18; }
		if($_[7] != -1){
			$$ref .= "$TEXT $black $tinyFont $x1 $y1 $_[7]\n";
		}
		if($_[8] != -1){
			$$ref .= "$TEXT $black $tinyFont $x2 $y1 $_[8]\n";
		}
	}
}

#*****************
sub find_in_list
#*****************
{
my $start = shift;
my $ref = shift;
my $N = shift;
my $pos = shift;

if(!$pos){ $pos=0; }
my $n = $N;
my $j = -1;
my $i = int($n/2);
while(1){
	my $start1=$ref->[$i]->[$pos];
	if($start1==$start){
		last;
	}
	elsif($start > $start1){
		if($n-$i<=1){ last; }
		$j = $i;
		$i = int(($n+$i)/2);
	}
	else{
		if($i-$j<=1){ --$i; last; }
		$n = $i;
		$i = int(($j+$n)/2);
	}
}
return $i;
}

#**********************************************
sub  normal_distribution
#**********************************************
{
my $x = shift;
my $a1=-1.26551223;
my $a2= 1.00002368;
my $a3= 0.37409196;
my $a4= 0.09678418;
my $a5=-0.18628806;
my $a6= 0.27886807;
my $a7=-1.13520398;
my $a8= 1.48851587;
my $a9=-0.82215223;
my $a10=0.17087277;
my $z = abs($x/sqrt(2));
my $t = 1.0 / (1.0 + 0.5 * $z);
my $y = $t*exp(-$z * $z + $a1 + $t * ($a2 + $t * ($a3 + $t * ($a4 + $t * ($a5 + $t *
     ($a6 + $t * ($a7 + $t * ($a8 + $t * ($a9 + $t * $a10)))))))));
if($x < 0.0){ $y = 2.0 - $y; }
$y = 1.0 - 0.5 * $y;
return $y;
}

#*****************************************
sub  parse_data
#*****************************************
{
my $data_ref = shift;
my $hashForm_ref = shift;

my @lines = split(/\r?\n\r?/,$$data_ref);
if(@lines==1){
	@lines = split(/\n?\r\n?/,$$data_ref);
}
my $iline = 0;
while($iline<@lines && $lines[$iline] =~ /^--/){
	my($name,$attribute,$filename) = ("","","");
	++$iline;
	while($iline<@lines && $lines[$iline]){
		my @items = split(/; /,$lines[$iline]);
		foreach my $item (@items){
			if($item !~ /=/){ next; }
			my($type,$value) = split(/=/,$item);
			if($type eq "name"){
				$value =~ s/"//g;
				$name = $value;
			}
			if($type eq "filename"){
				$value =~ s/"//g;
				$filename = $value;
			}
		}
		++$iline;
	}
	while($iline<@lines && !$lines[$iline]){ ++$iline; }
	while($iline<@lines && $lines[$iline] !~ /^--/){
		$attribute .= "$lines[$iline]\n";
		++$iline;
	}
	$attribute =~ s/\n+$//;
	if($name && $attribute){
		$name = substitute_chars($name);
		if($filename){
			$filename = substitute_chars($filename);
			$hashForm_ref->{$name} = [$filename,$attribute];
		}else{
			$attribute = substitute_chars($attribute);
			$hashForm_ref->{$name} = $attribute;
		}
	}
}
return;
}

#**************************************
sub  get_outputID
#**************************************
{
my $block = shift;

if(!$block){ $block=1; }
my @letters = ('1'..'9');
my $outputID=$letters[rand@letters];
push(@letters,'0');
for(my $i=0; $i<14; $i++){
	$outputID .= $letters[rand@letters];
}
for(my $id=$outputID; $id<$outputID+$block; ++$id){
	my $list = `ls -l $PATH_OUTPUT/$id.*`;
	if($list){
		`rm $PATH_OUTPUT/$id.*`;
	}
}
return $outputID;
}

#**************************************
sub   pattern_to_PFM
#**************************************
{
my $line_ref = shift;

my %hashCode=();
for(my $i=1; $i<@codes; ++$i){
	my @num=(0,0,0,0);
	my $ii = $i;
	my $sum = 0;
	for(my $k=0; $k<4; ++$k){
		if($ii%2){
			$num[$k]=1;
			$sum += 1;
		}
		$ii = int($ii/2);
	}
	for(my $k=0; $k<4; ++$k){
		$num[$k] = int(100*$num[$k]/$sum);
	}
	$hashCode{$codes[$i]} = join("\t",@num);
}
my $count;
for(my $i=0; $i<@$line_ref; ++$i){
	$line_ref->[$i] =~ s/^>[ \t]+/>/;
	my $line = $line_ref->[$i];
	if($line =~ /^>/){
		my($name,$pattern) = split(/\t/,$line);
		if(!$pattern){ print "No pattern found in line $i<br>"; return(0); }
		$pattern = uc($pattern);
		my @ch = split(//,$pattern);
		for(my $j=0; $j<@ch; ++$j){
			my $freq = $hashCode{$ch[$j]};
			if(!$freq){ print "Unrecognized pattern $pattern in line $i<br>"; return(0); }
			$line_ref->[$i] .= "\n$j\t$freq";
		}
		$line_ref->[$i] .= "\n";
		++$count
	}elsif($line !~ /^Headers/i && $line !~ /^Parameters/i){
		splice(@$line_ref,$i,1);
		--$i;
	}
}
if(!$count){ print "No patterns<br>"; return(0); }
return(1);
}

#**************************************
sub   check_upload
#**************************************
{
my $line_ref = shift;
my $file_type = shift;

my $count=0;
my $junk;
my $pattern_insert = 0;
my @headers=();
my %hashNames;
for(my $iline=0; $iline<@$line_ref; ++$iline){
	$line_ref->[$iline] =~ s/^>[ \t]+/>/;
	my $line = $line_ref->[$iline];
	while(!$line && $file_type !~ /motif/ && $iline<@$line_ref){
		splice(@$line_ref,$iline,1);
		$line_ref->[$iline] =~ s/^>[ \t]+/>/;
		$line = $line_ref->[$iline];
	}
	if($line =~ /^Headers|^Parameters/){
		if($line =~ /^Headers/ && $file_type eq "motif"){
			($junk,@headers)=split(/\t/,$line);
			if($headers[1] ne "Pattern"){
				splice(@headers,1,0,"Pattern");
				$pattern_insert = 1;
				$line_ref->[$iline] = join("\t","Headers:",@headers);
			}
		}
		next;
	}
	if($file_type eq "fasta"){
		if($line =~ /^>/){
			my($name,@list) = split(/\t/,$line);
			if($name eq ">"){
				error_message("Error in FASTA file: No sequence name in line $iline");
			}
			if($name =~ / /){
				$name =~ s/ /_/g;
				$line_ref->[$iline] = $name;
				if(@list){ $line_ref->[$iline] .= "\t".join("\t",@list); }
			}
			if($hashNames{$name}){
				error_message("Error in FASTA file: Duplicated sequence name in line $iline");
			}
			$hashNames{$name} = 1;
			++$count;
		}else{
			$line_ref->[$iline] =~ s/\s+//g;
			$line_ref->[$iline] =~ s/U/T/g;
			$line_ref->[$iline] =~ s/u/t/g;
			if($iline<2000 && length($line_ref->[$iline])>0 && $line_ref->[$iline] !~ /^[ACGTNacgtn]+$/){
				error_message("Error in FASTA file: wrong nucleotide in line $iline");
			}
		}
	}
	if($file_type eq "repeat" && $iline<1000){
		if($line =~ /^>/){ error_message("Error in line $iline"); }
		my($gap,$ind,$ratio,$len) = split(/\t/,$line);
		if($gap<0 || $gap>20 || $ratio<1 || !$len){ error_message("Error in line $iline"); }
		$line = $line_ref->[++$iline];
		my @items = split(/\s+/,$line);
		if(@items != $len*4){ error_message("Error in REPEAT file in line $iline"); }
		++$count;
	}
	if($file_type eq "search" && $iline<1000){
		if($line =~ /^>/){ error_message("Error in line $iline"); }
		my($motifName,$seqName,$strand,$len,$pos,$score) = split(/\t/,$line);
		++$count;
		if($motifName eq "NONE"){ next; }
		if(!$score || !$len || $strand !~ /^[\+-]$/){ error_message("Error in SEARCH file line $iline"); }
	}
	if($file_type eq "motif"){
		if(!@headers){
			$pattern_insert = 1;
			splice(@$line_ref,0,0,"Headers:\tName\tPattern");
			@headers = ("Name","Pattern");
			++$iline;
		}
		while($line !~ /^>/){
			splice(@$line_ref,$iline,1);
			if($iline==@$line_ref){ last; }
			$line = $line_ref->[$iline];
		}
		if($iline==@$line_ref){ last; }
		my @first_line_items = split(/\t/,$line);
		my $i1 = $iline;
		my $npos = 0;
		my $consensus="";
		while(1){
			++$i1;
			if($i1>=@$line_ref){
				$line_ref->[--$i1] .= "\n\n";
				last;
			}
			if($line_ref->[$i1] =~ /^>/){
				$line_ref->[--$i1] .= "\n";
				last;
			}
			$line_ref->[$i1] =~ s/^\s+//;
			$line = $line_ref->[$i1];
			if(length($line) < 2){ last; }
			if($line =~ /[A-z]/){ next; }
			my @m = split(/\t/,$line);
			my $ind = $npos;
			if(@m > 4){
				shift(@m);
			}
			my $sum = 0;
			for(my $i2=0; $i2<4; ++$i2){
				if($m[$i2] < 0){  error_message("Negative matrix line $i1"); }
				$sum += $m[$i2];
			}
			if(!$sum){  error_message("No data found in matrix line $i1"); }
			my $max=0;
			my $coef=1;
			my $x=0;
			for(my $i2=0; $i2<4; ++$i2){
				if($sum<=1.01){	$m[$i2] = int(100*$m[$i2]/$sum); }
				if($max<$m[$i2]){ $max=$m[$i2]; }
			}
			for(my $i2=0; $i2<4; ++$i2){
				my $y = 1;
				if($m[$i2] < 0.5*$max){ $y=0; }
				$x += $y*$coef;
				$coef *= 2;
			}
			$consensus .= $codes[$x];
			$line_ref->[$i1] = join("\t",$ind,@m);
			++$npos;
		}
		if($pattern_insert){
			splice(@first_line_items,1,0,$consensus);
			$line_ref->[$iline] = join("\t",@first_line_items);
		}
		$iline = $i1;
		++$count;
	}
	if($file_type eq "coordinates" && $iline<500){
		if($iline==0){ next; }
		if($line =~ /^>/){ error_message("Error in line $iline"); }
		my($Rname,$chr,$strand,$len,$pos) = split(/\t/,$line);
		++$count;
		if(!$pos || !$len || $strand !~ /^[\+-]$/ || $chr !~ /chr/){
			my ($chr,$start,$end) = split(/\t/,$line);
			if(!$start || $end<=$start || $chr !~ /chr/){
				error_message("Error in line $iline");
			}
		}
	}
	if($file_type eq "attributes" && $iline<1000){
		++$count;
	}
	if($file_type eq "conservation" && $iline<1000){
		if($line =~ /^>/){ ++$count; }
		elsif($line =~ /[A-Za-z;\t]/){ error_message("Error in conservation file, line $iline"); }
	}
	if($file_type eq "subset" && $iline<1000){
		if($line =~ /^>/){ error_message("Error in SUBSET file line $iline"); }
		++$count;
	}
}
if(!$count){ error_message("No records in $file_type file"); }
return(1);
}

#**************************************
sub  plot_frequency_distribution
#**************************************
{
my $inputFile = shift;
my $motif = shift;
my $outputID = shift;

open (INFO, $inputFile) or error_message("Input file not found.");
my $line = <INFO>;
$line =~ s/\n$//;
my ($junk,$junk,@headers)=split(/\t/,$line);
while(my $line = <INFO>){
	$line =~ s/\n$//;
	my ($name,$type,@data)=split(/\t/, $line);
	if($name eq $motif){
		my @dataAll = (\@data);
		$line = <INFO>;
		$line =~ s/\n$//;
		my ($name1,$type1,@data1)=split(/\t/, $line);
		push(@dataAll,\@data1);
		plot_data($name,\@headers,\@dataAll,$outputID+1);
		for(my $i=0; $i<@data; ++$i){
			$data[$i] += $data1[$i];
		}
		pop(@dataAll);
		for(my $i=0; $i<3; ++$i){
			$line = <INFO>;
			$line =~ s/\n$//;
			my ($name2,$type2,@data2)=split(/\t/, $line);
			push(@dataAll,\@data2);
		}
		plot_data($name,\@headers,\@dataAll,$outputID);
		last;
	}
}
return;
}

#*********************************
sub  plot_data
#*********************************
{
my $name = shift;
my $x_ref = shift;
my $data_ref = shift;
my $outputID = shift;

my @colors = ($magenta,$blue,$green,$gray);

# Find minimum and maximum
my $minx=$x_ref->[0];
my $maxx=$x_ref->[@$x_ref-1];
my $miny = 1000000;
my $maxy = -1000000;
foreach my $y_ref (@$data_ref){
	foreach my $y (@$y_ref){
		if($miny > $y) { $miny=$y; }
		if($maxy < $y) { $maxy=$y; }
	}
}
if($miny > 0){ $miny=0; }
if($miny > $maxy){ $miny=0; $maxy=1; }
my $spanx = $maxx-$minx; if($spanx <= 0) { $spanx=1; }
my $spany = $maxy-$miny; if($spany <= 0) { $spany=0.1; }
my $deltax = 20;
if($spanx > 200) { $deltax = 50; }
if($spanx > 500) { $deltax = 100; }
if($spanx > 1000) { $deltax = 200; }
if($spanx > 2000) { $deltax = 500; }
if($spanx > 5000) { $deltax = 1000; }
if($spanx > 10000) { $deltax = 2000; }
if($spanx > 20000) { $deltax = 5000; }
if($spanx > 50000) { $deltax = 10000; }
my $deltay = 0.0001;
if($spany > 0.001) { $deltay = 0.0002; }
if($spany > 0.002) { $deltay = 0.0005; }
if($spany > 0.005) { $deltay = 0.001; }
if($spany > 0.01) { $deltay = 0.002; }
if($spany > 0.02) { $deltay = 0.005; }
if($spany > 0.05) { $deltay = 0.01; }
if($spany > 0.1) { $deltay = 0.02; }
if($spany > 0.2) { $deltay = 0.05; }
if($spany > 0.5) { $deltay = 0.1; }
if($spany > 1) { $deltay = 0.2; }
if($spany > 2) { $deltay = 0.5; }
if($spany > 5) { $deltay = 1; }
if($spany > 10) { $deltay = 2; }
if($spany > 20) { $deltay = 5; }
if($spany > 50) { $deltay = 10; }
if($spany > 100) { $deltay = 20; }
if($spany > 200) { $deltay = 50; }
if($spany > 500) { $deltay = 100; }
if($spany > 1000) { $deltay = 200; }
if($spany > 2000) { $deltay = 500; }
if($spany > 5000) { $deltay = 1000; }

open (OUT1, ">$PATH_OUTPUT/$outputID.txt") or error_message("Output file not opened.");
my $WID = 400;
my $HGT = 150;
my $scalex = 0.875*$WID/$spanx;
my $scaley = 0.767*$HGT/$spany;
my $x_center = int(0.400*$WID);
my $y_center = int(0.333*$HGT);
my $nobjects = 1000+@$x_ref*10;
print OUT1 "1\n$x_center\n$y_center\n$nobjects\n";

# Axis horiz & vertical
my $border = 6;
my $x1 = -$border;
my $x2 = int(($maxx-$minx)*$scalex) + $border;
my $y1 = 0;
my $y2 = int(($maxy-$miny)*$scaley);
print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x2 $y1\n";
print OUT1 "$LINE $black $thin $solid  $x1 $y1 $x1 $y2\n";

# Make ticks on the horizontal axis
my $tick = 5;
if($HGT <= 80 || $WID <= 80){ $tick = 2; }
my $y3 = $y1-$tick;
my $y4 = $y1-15;
for(my $i = 0; $i <= $maxx/$deltax; ++$i){
	my $x = int($i*$deltax);
	if($x < $minx){ next; }
	my $x3 = int(($x-$minx)*$scalex);
	print OUT1 "$LINE $black $thin $solid  $x3 $y1 $x3 $y3\n";
	if($HGT > 80 && $WID > 80){
		print OUT1 "$TEXT $black $smallFont $x3 $y4 $x\n";
	}
}
for(my $i = 1; $i < -$minx/$deltax; ++$i) {
	my $x = -int($i*$deltax);
	my $x3 = int(($x-$minx)*$scalex);
	print OUT1 "$LINE $black $thin $solid  $x3 $y1 $x3 $y3\n";
	if($HGT > 80 && $WID > 80){
		print OUT1 "$TEXT $black $smallFont $x3 $y4 $x\n";
	}
}
my $x5 = (($maxx + $minx)/2-$minx)*$scalex;
my $y5 = -30;

# Make ticks on the vertical axis
my $x3 = $x1-$tick;
my $x4 = $x1-20;
for(my $i = 0; $i <= $maxy/$deltay; ++$i) {
	my $y = int($i*$deltay);
	if($y < $miny){ next; }
	my $y3 = int(($y-$miny)*$scaley);
	my $y4 = $y3+2;
	print OUT1 "$LINE $black $thin $solid  $x1 $y3 $x3 $y3\n";
	if($HGT > 80 && $WID > 80){
		print OUT1 "$TEXT $black $smallFont $x4 $y4 $y\n";
	}
}
for (my $i = 1; $i < -$miny/$deltay; ++$i) {
	my $y = -int($i*$deltay);
	my $y3 = int(($y-$miny)*$scaley);
	my $y4 = $y3+4;
	print OUT1 "$LINE $black $thin $solid  $x1 $y3 $x3 $y3\n";
	if($HGT > 80 && $WID > 80){
		print OUT1 "$TEXT $black $smallFont $x4 $y4 $y\n";
	}
}
my $x5 = 20;
my $y5 = int(($maxy-$miny)*$scaley) + 20;

# Print lines
foreach (my $if=0; $if<@$data_ref; ++$if){
	my $y_ref = $data_ref->[$if];
	for(my $i=1; $i<@$x_ref; ++$i){
		my $x0 = int(($x_ref->[$i-1]-$minx)*$scalex);
		my $y0 = int(($y_ref->[$i-1]-$miny)*$scaley);
		my $x1 = int(($x_ref->[$i]-$minx)*$scalex);
		my $y1 = int(($y_ref->[$i]-$miny)*$scaley);
		print OUT1 "$LINE $colors[$if] 1 $solid  $x0 $y0 $x1 $y1\n";
		#print "$i $x_ref->[$i] $y_ref->[$i]<br>\n";
	}
}
close (OUT1);
`./togif.exe $PATH_OUTPUT/$outputID.txt $PATH_OUTPUT/$outputID.gif raster.par palette.par $WID $HGT`;
return;
}

#*****************************************
sub  get_array_lists
#*****************************************
{
my $ref = shift;

my $items="";
my $descriptions="";
my $files="";
foreach my $ref1 (@$ref){
	my($item,$descr,$file) = @$ref1;
	if(!$items){ $items = "\"".$item."\""; }
	else{ $items .= ",\"".$item."\""; }
	if(!$descriptions){ $descriptions = "\"".$descr."\""; }
	else{ $descriptions .= ",\"".$descr."\""; }
	if(!$files){ $files = "\"".$file."\""; }
	else{ $files .= ",\"".$file."\""; }
}
return($items,$descriptions,$files);
}

#*****************************************
sub   file_split
#*****************************************
{
my $file_split = shift;
my $file_split_result = shift;
my $file_split_description = shift;
my $lines_ref = shift;

if(!$lines_ref){ error_message("Empty region coordinates!"); }
my $file_split_regions = $lines_ref->[0];
if(@$lines_ref>1){ $file_split_regions = join(';',@$lines_ref); }
$file_split_regions =~ s/\t|\n/;/g;
$file_split_description =~ s/\t|\n/ /g;
my @items = split(/[;, ]+/,$file_split_regions);
my @regions=();
foreach my $item (@items){
	my($start,$end) = split(/-/,$item);
	if($start >= $end || $start<0 || $end<=0){ error_message("Wrong region coordinates!"); }
	push(@regions,[$start,$end]);
}
my $destinFile = "$PATH_DATA/$loginname-$file_split_result";
my $sourceFile = "$PATH_DATA/$loginname-$file_split";
if($file_split =~ /^public-/){
	$sourceFile = "$PATH_DATA/$file_split";
}
my $Tname;
my $sequence;
if(!open (OUTPUT, ">$destinFile")){ error_message("Cannot write to file"); }
if(!open (INFO, $sourceFile)){ error_message("Source file not found!"); }
my $count = 0;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^>/){
		if($Tname && $sequence){
			$count += split_sequence($Tname,$sequence,\@regions);
		}
		my @items = split(/\s+/,$line);
		$Tname = $items[0];
		$Tname =~ s/^>//;
		$sequence = "";
	}else{
		$line =~ s/\s//g;
		if($line){
			$sequence .= $line;
		}
	}
}
if ($Tname && $sequence){
	$count += split_sequence($Tname,$sequence,\@regions);
}
close INFO;
close OUTPUT;
if(!$count){ error_message("Output file is empty; apparently coordinates are out of bounds."); }

#Split conservation file
my $destinFile1 = $destinFile;
my $sourceFile1 = $sourceFile;
$destinFile1 =~ s/\.fa$/.cons/;
$sourceFile1 =~ s/\.fa$/.cons/;
my $Tname;
my @conserv=();
my $interval=1;
if(open (INFO, $sourceFile1)){
	if(!open (OUTPUT, ">$destinFile1")){ error_message("Cannot write to file"); }
	while(my $line = <INFO>){
		if($line =~ /^Parameters|^Headers/){
			print OUTPUT $line;
			if($line =~ /^Parameters/){
				my %hash;
				read_config_line($line,\%hash);
				$interval = $hash{"interval"};
				if($interval<1){ $interval=1; }
			}
			next;
		}
		$line =~ s/\n$//;
		if($line =~ /^>/){
			if ($Tname && $sequence){
				split_conservation($Tname,\@conserv,\@regions,$interval);
			}
			my @items = split(/\t/,$line);
			$Tname = $items[0];
			$Tname =~ s/^>//;
			@conserv=();
		}else{
			$line =~ s/\s//g;
			if($line){
				push(@conserv,split(/,/,$line));
			}
		}
	}
	if ($Tname && $sequence){
		split_conservation($Tname,\@conserv,\@regions,$interval);
	}
	close INFO;
	close OUTPUT;
}

#Generate coordinates file
my $destinFile1 = $destinFile;
my $sourceFile1 = $sourceFile;
$destinFile1 =~ s/\.fa$/.coord/;
$sourceFile1 =~ s/\.fa$/.coord/;
if(open (INFO, $sourceFile1)){
	if(!open (OUTPUT, ">$destinFile1")){ error_message("Cannot write to file"); }
	while(my $line = <INFO>){
		if($line =~ /^Parameters|^Headers/){
			print OUTPUT $line;
			next;
		}
		$line =~ s/\n$//;
		if(!$line){ next; }
		my($Tname,$chr,$strand,$size,$start) = split(/\t/,$line);
		foreach my $ref (@regions){
			my($start1,$end1) = @$ref;
			my $newname = $Tname;
			if(@regions>1){
				$newname = "$Tname-$start1-$end1";
			}
			if($end1 > $size){ $end1 = $size; }
			if($end1 <= $start1){ next; }
			my $newsize = $end1 - $start1;
			my $newstart = $start + $start1;
			if($strand eq "-"){ $newstart = $start + $size - $end1; }
			print OUTPUT "$newname\t$chr\t$strand\t$newsize\t$newstart\n";
		}
	}
	close INFO;
	close OUTPUT;
}
#Generate attribute file
my $destinFile1 = $destinFile;
my $sourceFile1 = $sourceFile;
$destinFile1 =~ s/\.fa$/.attr/;
$sourceFile1 =~ s/\.fa$/.attr/;
if(open (INFO, $sourceFile1)){
	if(!open (OUTPUT, ">$destinFile1")){ error_message("Cannot write to file"); }
	while(my $line = <INFO>){
		if($line =~ /^Parameters|^Headers/){
			print OUTPUT $line;
			next;
		}
		$line =~ s/\n$//;
		if(!$line){ next; }
		my($Tname,@items) = split(/\t/,$line);
		foreach my $ref (@regions){
			my($start1,$end1) = @$ref;
			my $newname = $Tname;
			if(@regions>1){
				$newname = "$Tname-$start1-$end1";
			}
			print OUTPUT "$newname\t".join("\t",@items)."\n";
		}
	}
	close INFO;
	close OUTPUT;
}
open(OUT, ">$PATH_INFO/$loginname-config1.txt");
open(INFO, "$PATH_INFO/$loginname-config.txt");
my %hash=();
while(my $line = <INFO>){
	if($line !~ /^type_fasta=$file_split_result/){ print OUT $line; }
	if($line =~ /^type_fasta=$file_split/){
		read_config_line($line,\%hash);
	}
}
close INFO;
print OUT "type_fasta=$file_split_result";
foreach my $key (%hash){
	my $newFilename = $file_split_result;
	if($key eq "coordinates"){ $newFilename =~ s/\.fa$/\.coord/; }
	elsif($key eq "attributes"){ $newFilename =~ s/\.fa$/\.attr/; }
	elsif($key eq "conservation"){ $newFilename =~ s/\.fa$/\.cons/; }
	else{ next; }
	print OUT "\t$key=$newFilename";
}
print OUT "\tdescription=$file_split_description\n";
close OUT;
system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
system("rm $PATH_INFO/$loginname-config1.txt");
return;
}

#*****************************************
sub   split_sequence
#*****************************************
{
my $Tname = shift;
my $sequence = shift;
my $region_ref = shift;

my $length = length($sequence);
my $count = 0;
foreach my $ref (@$region_ref){
	my($start,$end) = @$ref;
	my $newname = $Tname;
	if(@$region_ref>1){
		$newname = "$Tname-$start-$end";
	}
	if($end > $length){ $end = $length; }
	if($end<=$start){ next; }
	++$count;
	my $subseq = substr($sequence,$start,$end-$start);
	print OUTPUT ">$newname\n";
	my $len = length($subseq);
	for(my $i=0; $i<$len; $i+=100){
		my $short = $len - $i;
		if($short > 100){ $short=100; }			
		my $string = substr($subseq,$i,$short);
		print OUTPUT "$string\n";
	}
}
return $count;
}

#*****************************************
sub   split_conservation
#*****************************************
{
my $Tname = shift;
my $cons_ref = shift;
my $region_ref = shift;
my $interval = shift;

my $length = @$cons_ref;
foreach my $ref (@$region_ref){
	my($start,$end) = @$ref;
	my $newname = $Tname;
	if(@$region_ref>1){
		$newname = "$Tname-$start-$end";
	}
	$start = int($start/$interval);
	$end = int($end/$interval);
	if($end > $length){ $end = $length; }
	if($end<=$start){ next; }
	print OUTPUT ">$newname\n";
	my $len = $end - $start;
	for(my $i=0; $i<$len; $i+=100){
		my $end1 = $start+$i+100;
		if($end1 > $length){ $end1=$length; }
		print OUTPUT $cons_ref->[$start+$i];
		for(my $j=$start+$i+1; $j<$end1; ++$j){
			print OUTPUT ",$cons_ref->[$j]";
		}
		print OUTPUT "\n";
	}
}
return;
}

#*****************************************
sub   advanced_search
#*****************************************
{
my $seq_ref = shift;
my $hashSeq_ref = shift;
my $header_line = shift;

my $search_term = $hashInput{"search_term"};
my $field = $hashInput{"field"};
my $distanceTSS = $hashInput{"distanceTSS"};
my $distance = $hashInput{"distance"};
my $combineMotifs = $hashInput{"combineMotifs"};

my @TFlist;
for(my $iTF=1; $iTF<=3 && $hashInput{"TF$iTF"} && $hashInput{"TF$iTF"} ne "none"; ++$iTF){
	my $TF = $hashInput{"TF$iTF"};
	if(!$TF){ last; }
	my $startTF=$hashInput{"startTF$iTF"};
	my $endTF=$hashInput{"endTF$iTF"};
	my $nhitsTF=$hashInput{"nhitsTF$iTF"};
	my $conservTF=$hashInput{"conservTF$iTF"};
	my $scoreTF=$hashInput{"scoreTF$iTF"};
	my $strandTF=$hashInput{"strandTF$iTF"};
	push(@TFlist,[$TF,$startTF,$endTF,$nhitsTF,$conservTF,$scoreTF,$strandTF]);
}
# First step: find sequences
my %hashSelect;
my ($no,$name,@headers)=split(/\t/,$header_line);
my ($isymbol,$idist) = (-1,-1);
for(my $i=0; $i<@headers; ++$i){
	if($headers[$i] =~ /symbol/i){ $isymbol=$i; }
	if($headers[$i] =~ /distance/i){ $idist=$i; }
}
if($search_term){
	foreach my $name (keys %$hashSeq_ref){
		my $ref = $hashSeq_ref->{$name};
		if($field =~ /default|sequence/){
			if($name eq $search_term){ $hashSelect{$name}=1; }
		}
		if($field =~ /default|symbol/ && $isymbol >=0 && $ref->[$isymbol] =~ /$search_term/){
			my @symbol =split(/[,;]/,$ref->[$isymbol]);
			for(my $i=0; $i<@symbol; ++$i){
				if(uc($symbol[$i]) eq uc($search_term)){
					if($distanceTSS && $idist>=0){
						my @dist =split(/[,;]/,$ref->[$idist]);
						if(abs($dist[$i]) <= $distanceTSS){ $hashSelect{$name}=1; }
					}else{
						$hashSelect{$name}=1;
					}
					last;
				}
			}
		}
	}
	if(!%hashSelect){ not_found(); }
}
elsif($distanceTSS && $idist>=0){
	foreach my $name (keys %$hashSeq_ref){
		my $ref = $hashSeq_ref->{$name};
		my @dist =split(/[,;]/,$ref->[$idist]);
		for(my $i=0; $i<@dist; ++$i){
			if(abs($dist[$i]) <= $distanceTSS){
				$hashSelect{$name}=1;
				last;
			}
		}
	}
	if(!%hashSelect){ not_found(); }
}
if(!@TFlist){
	if(!%hashSelect){
		return;
	}
	@$seq_ref = sort keys %hashSelect;
	return;
}
# Second step: find motifs
my $file_search = $hashInput{"file_search"};
my $file_search_output = $hashInput{"file_search_output"};
my $fileSearch = "$PATH_DATA/$file_search";
if($file_search !~ /^public-/ && !$hashInput{"data_saved"}){ $fileSearch = "$PATH_DATA/$loginname-$file_search"; }
if($file_search_output){ $fileSearch = "$PATH_OUTPUT/$file_search_output"; }
if(!open(INFO, $fileSearch)){
	error_message("Search output file not found!");
}
my %hashTF=();
for(my $i=0; $i<@TFlist; ++$i){
	$hashTF{$TFlist[$i]->[0]} = $i;
}
my %hashTFnum;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line =~ /^Headers|^Parameters/){ next; }
	my($motifName,$seqName,$strand,$len,$pos,$score,$pattern,$cons) = split(/\t/,$line);
	if(%hashSelect && !$hashSelect{$seqName}){ next; }
	my $iTF = $hashTF{$motifName};
	if(!defined($iTF)){ next; }
	if($pos<$TFlist[$iTF]->[1] || $pos>=$TFlist[$iTF]->[2] || $cons ne "" && $cons<$TFlist[$iTF]->[4] || $score<$TFlist[$iTF]->[5]){ next; }
	if($strand eq "+"){ $strand=1; }
	else{ $strand=-1; }
	if($strand*$TFlist[$iTF]->[5] < 0){ next; }
	push(@{$hashTFnum{$seqName}->[$iTF]},$pos);
}
@$seq_ref=();
my $combine = $hashInput{"combineMotifs"};
my $distOption = $hashInput{"distance"};
foreach my $seqName (sort keys %hashTFnum){
	my $ref = $hashTFnum{$seqName};
	if(!$ref){ next; }
	my $selectedAnd=1;
	my $selectedOr=0;
	my @positionList=();
	for(my $iTF=0; $iTF<@TFlist; ++$iTF){
		my $ref1 = $ref->[$iTF];
		if(!$ref1){ $selectedAnd=0; next; }
		my $n1 = @$ref1;
		if($n1 < $TFlist[$iTF]->[3]){ $selectedAnd=0; next; }
		$selectedOr=1;
	}
	if($combine =~ /and|dist/ && !$selectedAnd || $combine eq "or" && !$selectedOr){ next; }
	if($combine ne "dist" || @TFlist==1){
		for(my $iTF=0; $iTF<@TFlist; ++$iTF){
			my $ref1 = $ref->[$iTF];
			if($ref1){ push(@positionList,@$ref1); }
		}
	}else{
		for(my $iTF=0; $iTF<@TFlist; ++$iTF){
			my $ref1 = $ref->[$iTF];
			my %hashHits=();
			foreach my $x (@$ref1){ $hashHits{$x}=1; }
			for(my $jTF=$iTF+1; $jTF<@TFlist; ++$jTF){
				my $ref2 = $ref->[$jTF];
				my %hashHits1=%hashHits;
				foreach my $x (@$ref2){ $hashHits1{$x} += 2; }
				my @sorted = (sort {$a<=>$b} keys %hashHits1);
				for(my $k=0; $k<@sorted; ++$k){
					my $x = $sorted[$k];
					for(my $l=$k+1; $l<@sorted && $sorted[$l]-$x<$distOption; ++$l){
						my $x1 = $sorted[$l];
						if($hashHits1{$x1} != $hashHits1{$x} || $hashHits1{$x}==3 && $hashHits1{$x1}==3){
							push(@positionList,int(($x+$x1)/2));
						}
					}
				}
			}
		}
	}
	if(@positionList){
		push(@$seq_ref,$seqName);
	}
}
@$seq_ref = sort @$seq_ref;
return;
}

#*****************************************
sub   subset_table
#*****************************************
{
my $source = shift;
my $filename = shift;
my $list_ref = shift;
if(!open(SOURCE,$source)){ error_message("File not found!"); }
if(!open(DESTIN,">$filename")){ error_message("Cannot write to file!"); }
my %hash=();
foreach my $item (@$list_ref){
	$item =~ s/\t.+//;
	$hash{$item}=1;
}
while(my $line = <SOURCE>){
	$line =~ s/\n$//;
	my @items = split(/\t/,$line);
	if($line =~ /^Parameters|^Headers/ || $hash{$items[0]}){
		print DESTIN "$line\n";
	}
}
close SOURCE;
close DESTIN;
return 1;
}

#*****************************************
sub   subset_fasta
#*****************************************
{
my $source = shift;
my $filename = shift;
my $list_ref = shift;
if(!open(SOURCE,$source)){ error_message("File $source not found!"); }
if(!open(DESTIN,">$filename")){ error_message("Cannot write to file $filename"); }
my %hash=();
foreach my $item (@$list_ref){
	$item =~ s/\t.+//;
	$hash{$item}=1;
}
my $include = 1;
my $count = 0;
while(my $line = <SOURCE>){
	if($line =~ /^>/){
		my $name = $line;
		$name =~ s/\n$//;
		$name =~ s/^>//;
		if($hash{$name}){ $include = 1; $count++; }
		else{ $include = 0; }
	}
	if($include){
		print DESTIN $line;
	}
}
close SOURCE;
close DESTIN;
if(!$count){ error_message("Output file empty; apparently sequence names do not match."); }
return 1;
}

#**************************************
sub  extract_sequence_from_genome
#**************************************
{
my $lines = shift;
my $filename = shift;
my $genome = shift;

my $length = 0;
my $bed_format = 0;
my ($chr,$start,$end) = split(/\t/,$lines->[1]);
if($chr =~ /^chr/ && $start>=0 && $end > $start){
	$bed_format = 1;
}
my $num = 1;
for(my $i=1; $i<@$lines; ++$i){
	my ($name,$chr,$strand,$size,$start,$end);
	if($bed_format){
		($chr,$start,$end,$name) = split(/\t/,$lines->[$i]);
		if(!$name){
			$name = sprintf("R%07d",$num);
			++$num;
		}
		$start--;
		$size = $end-$start;
	}else{
		($name,$chr,$strand,$size,$start) = split(/\t/,$lines->[$i]);
	}
	if($chr !~ /^chr/ || $strand !~ /-|\+/ || $size<=0 || $start<0){
		error_message("ERROR in genome coordinates, line $i.<p>$lines->[$i]\n");
	}
	$length += $size;
}
if($length > 100000000){
	error_message("ERROR: Total sequence length ($length) is more than 100 Kb");
}
if($filename =~ /\./){ $filename =~ s/\..*$/.coord/; }
else{ $filename .= ".coord"; }
my $destinFile = "$PATH_DATA/$loginname-$filename";
if(!open(DESTIN,">$destinFile")){ error_message("Cannot write to file $filename"); }
foreach my $line (@$lines){
	print DESTIN "$line\n";
}
close DESTIN;
print "<HTML><HEAD><TITLE>Extracting sequence</TITLE>\n";
print_header();
print "<h3>Sequences are being extracted from the genome</h3>\n";
print "This will take 10-20 min. After the process is finished, the file will appear in the\n";
print "list of sequence files. However, you will need to refresh the screen of the main menu\n";
print "to see it. There will be no notification that the process is finished.\n";
print "Now you can press the 'Continue' button to return to the main menu.<p>\n";
print "<FORM NAME=cisfinder ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
print "<INPUT NAME=\"loginname\" TYPE=\"hidden\" VALUE=\"$loginname\">\n";
print "<INPUT NAME=\"sessionID\" TYPE=\"hidden\" VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=\"continue\">\n";
print "<INPUT NAME=\"continue\" TYPE=\"submit\" VALUE=Continue>\n";
print "</FORM>\n";
print "</BODY>\n";
print "</HTML>\n";
my $fastaFile = $filename;
$fastaFile =~ s/\.coord/.fa/;
my $configLine = "type_extracted=$fastaFile\tcoordinates=$filename";
my $descr = $hashInput{"description"};
if($descr){ $configLine .= "\tdescription=$descr"; }
open(OUT, ">>$PATH_INFO/$loginname-config.txt");
print OUT "$configLine\n";
close OUT;
my $pid = fork();
if($pid){
	exit(0);
}
my %hashSeq;
my %hashGenParts;
my %hashChr;
my $input_file = "$PATH_DATA/$loginname-$filename";
my $fasta_file = $filename;
$fasta_file =~ s/\.coord/.fa/;
my $output_file = "$PATH_DATA/$loginname-$fasta_file";
my $genome_path = "$PATH_PROG/genome/$genome";
if($genome ne "mm9" && $genome ne "hg19" && $genome ne "hg18"){
	print "ERR: wrong genome!\n";
	exit(0);
}
#Read input file
open (INFO, $input_file) or die $!;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	my($seqName,$chr,$strand,$blockSizes,$tStarts) = split(/\t/, $line);
		if(!$tStarts || !$blockSizes){ next; }
	$hashSeq{$seqName} = "$chr\t$strand\t$blockSizes\t$tStarts";
	my @tStart = split(/,/, $tStarts);
	my $Tstart = $tStart[0];
	my $part = int($Tstart/$CHROM_LENGTH);
	if(!$hashChr{$chr}){ $hashChr{$chr} = 1; }
	if($hashGenParts{$chr}<$part+1){ $hashGenParts{$chr} = $part+1; }
}
close INFO;
# Loop for chromosome number
open (OUTPUT, ">$output_file") or die $!;
foreach my $chr_name (sort keys %hashChr){
	my $chr_filename = "$genome_path/$chr_name.fa";
	#if($chr_name ne "chr8"){ next; }
	print "$chr_name\n";
	if(!open (INFO, $chr_filename)){ print "not found!\n"; next; }
	my $sequence = "";
	for(my $ipart=0; $ipart<$hashGenParts{$chr_name}; ++$ipart){
		my $length = length($sequence);
		if($length > $CHROM_LENGTH){
			$sequence = substr($sequence,$CHROM_LENGTH,$length-$CHROM_LENGTH);
		}else{
			$sequence = "";
		}
		$length = length($sequence);
		my $offset = $CHROM_LENGTH*$ipart;
		my $seqEnd = $CHROM_LENGTH*($ipart+1);
		my $count = 0;
		while (my $line = <INFO>){
			$line =~ s/\n$//;
			if($line =~ />/){ next; }
			$line =~ s/\s//g;
			my $len = length($line);
			$length += $len;
			$sequence .= $line;
			if($length >= $CHROM_LENGTH+$BUFFER){ last; }
		}
		if(!$sequence){ last; }
		#print "$ipart $hashGenParts{$chr_name} $length\n";
		foreach my $seqName (keys %hashSeq){
			my ($chr,$strand,$blockSizes,$tStarts) = split(/\t/, $hashSeq{$seqName});
			if($chr ne $chr_name){ next; }
			if($tStarts < $offset || $tStarts >= $seqEnd){ next; }
			my $output = substr($sequence,$tStarts-$offset,$blockSizes);
			if($strand eq "-"){
				$output = reverse($output);
				$output =~ tr/ATGCatgc/TACGtacg/;
			}
			print_fasta($seqName,$output);
		}
	}
	close INFO;
}
close OUTPUT;
open(OUT, ">$PATH_INFO/$loginname-config1.txt");
open(INFO, "$PATH_INFO/$loginname-config.txt");
my $configLine;
while(my $line = <INFO>){
	if($line =~ /^type_extracted=$fasta_file/){
		$configLine = $line;
	}elsif($line !~ /^type_fasta=$fasta_file/){
		print OUT $line;
	}
}
$configLine =~ s/_extracted/_fasta/;
print OUT $configLine;
close INFO;
close OUT;
if($configLine){
	system("cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt");
}
system("rm $PATH_INFO/$loginname-config1.txt");
exit(0);
}

#*********************************
sub  print_fasta
#*********************************
{
my $header = shift;
my $sequence = shift;

print OUTPUT ">$header\n";
my $len = length($sequence);
for(my $i=0; $i<$len; $i+=100){
	my $short = $len - $i;
	if($short > 100){ $short=100; }			
	my $string = substr($sequence,$i,$short);
	print OUTPUT "$string\n";
}
return;
}

#**************************************
sub  error_message
#**************************************
{
my $message = shift;
my $comment = shift;

my ($header,@lines) = split(/\n/,$message);
my $text = "<H3>ERROR: $header</H3>\n";
foreach my $line (@lines){
	$text .= "$line<br>\n";
}
if($comment =~ /^\d+$/){
	`echo "ERROR: $message" >> $PATH_OUTPUT/$comment.txt`;
	exit(0);
}elsif($comment eq "register"){
	print "<HTML><HEAD><TITLE>CisFinder error message</TITLE>\n";
	print_header();
	print $text."<p>\n";
	print "&nbsp; &nbsp;  &nbsp; &nbsp; &nbsp; <A HREF=../index.html><IMG SRC=../images/button_home.gif BORDER=0></A><p>";
}else{
	terminal_window($text,$comment);
}
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub  get_header
#**************************************
{
my $onload = shift;

my $text = "";
$text .= "</HEAD><BODY BGCOLOR=white";
if($onload){
	$text .= " onLoad=\"$onload\"";
}
$text .= ">\n";
$text .= "<TABLE BGCOLOR=black>\n";
$text .= "<TR><TD WIDTH=168 background=../images/head3.gif><IMG SRC=../images/head1.gif BORDER=0 useMap=#mapHead1></A></TD>\n";
$text .= "<TD WIDTH=678 background=../images/head4.gif ALIGN=CENTER><IMG SRC=../images/head2.gif BORDER=0 ALT=\"CisFinder: DNA motif finder\" useMap=#mapHead2></TD></TR>\n";
$text .= "<map NAME='mapHead1'>\n";
$text .= "<area shape=rect href=https://www.nih.gov coords=0,0,82,50>\n";
$text .= "<area shape=rect href=https://www.irp.nia.nih.gov coords=83,0,167,50>\n";
$text .= "<area shape=rect href=https://www.irp.nia.nih.gov coords=0,51,167,86>\n";
$text .= "</map>\n";
$text .= "<map NAME='mapHead2'>\n";
$text .= "<area shape=rect href=../index.html coords=0,0,677,80>\n";
$text .= "</map>\n";
$text .= "</TABLE>\n";
return $text;
}

#**************************************
sub  print_header
#**************************************
{
my $onload = shift;
my $text = get_header($onload);
print $text;
return;
}

#**************************************
sub change_password_form
#**************************************
{
my $loginname = shift;
my $sessionID = shift;

print "<HEAD><title>Password Change Form</title>\n";
print "<SCRIPT LANGUAGE=\"JavaScript\"><!--\n";
print "function check_onsubmit() {\n";
print "	if(!document.passwd_change.passwd.value){\n";
print "		alert(\"You need to enter existing password\");\n";
print "		return false;\n";
print "	}\n";
print "	if(!document.passwd_change.new_passwd.value){\n";
print "		alert(\"You need to enter new password\");\n";
print "		return false;\n";
print "	}\n";
print "	if(document.passwd_change.new_passwd.value != document.passwd_change.new_passwd1.value){\n";
print "		alert(\"New passwords do not match. Please retype them.\");\n";
print "		return false;\n";
print "	}\n";
print "	document.passwd_change.submit();\n";
print "}\n";
print "// -->\n";
print "</SCRIPT></HEAD>\n";
print "<BODY BGCOLOR=WHITE>\n";
print "<H2>Change password</H2>\n";
print "<FORM NAME=passwd_change METHOD=POST ACTION=$CGI_ADDRESS/cisfinder.cgi>\n";
print "<TABLE BORDER=0>\n";
print "<TR><TD>Current password:<TD><input  type=\"password\" name=\"passwd\" value=\"\" style=width:200px;>\n";
print "<TR><TD>New password:<TD><input type=\"password\" name=\"new_passwd\" value=\"\" style=width:200px;>\n";
print "<TR><TD>Retype new password:<TD><input type=\"password\" name=\"new_passwd1\" value=\"\" style=width:200px;>\n";
print "<TR><TD><TD><INPUT TYPE=BUTTON VALUE = \"Change password\" onClick=\"check_onsubmit();\" style=width:200px;><p>\n";
print "<TR><TD><TD><INPUT TYPE=BUTTON VALUE = \"Cancel (close)\" onClick=\"window.close();\" style=width:200px;><p>\n";
print "</TABLE>\n";
print "<INPUT TYPE=hidden NAME=loginname VALUE=\"$loginname\">\n";
print "<INPUT TYPE=hidden NAME=sessionID VALUE=\"$sessionID\">\n";
print "<INPUT NAME=\"action\" TYPE=\"hidden\" VALUE=update_password>\n";
print "</FORM><p>\n";
exit(0);
}

#**************************************
sub update_password
#**************************************
{
my $loginname = shift;
my $new_passwd = $hashInput{"new_passwd"};

my @letters = ('A' .. 'Z', 'a' .. 'z', '0' .. '9', '_', '.');
open (INFO, "$PATH_INFO/login.txt") or die $!;
open (OUT, ">$PATH_INFO/login_temp.txt") or die $!;
while(my $line = <INFO>){
	chomp $line;
	my($name,$passwd,$lastname,$firstname,$email)=split(/\t/, $line);
	if($name eq $loginname){
		my $salt = $letters[rand@letters] . $letters[rand@letters];
 		my $encryptedPsw = $new_passwd;
		for(my $i=0; $i<47; $i++){
			$encryptedPsw = crypt $encryptedPsw, $salt;
		}
		$line = join("\t",$name,$encryptedPsw,$lastname,$firstname,$email);
	}
	print OUT "$line\n";
}
close INFO;
close OUT;
`cp $PATH_INFO/login_temp.txt $PATH_INFO/login.txt`;
`rm $PATH_INFO/login_temp.txt`;
terminal_window("<H2>Your password has been changed</H2>");
}

#**************************************
sub start_session
#**************************************
{
my $loginname = shift;
my $sessionID;
my @letters = ('A' .. 'Z', 'a' .. 'z', '0' .. '9', '_', '.');
my $n=@letters;
for(my $i=0; $i<12; $i++){
	$sessionID .= $letters[rand $n];
}
open(OUT, ">$PATH_INFO/$loginname-config1.txt");
open(INFO, "$PATH_INFO/$loginname-config.txt");
my $nLines;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if(!$line){ next; }
	$nLines++;
	if($line !~ /^sessionID=/){
		print OUT $line."\n";
	}
}
close INFO;
my $salt = $letters[rand@letters] . $letters[rand@letters];
my $encryptedID = $sessionID;
for(my $i=0; $i<47; $i++){
	$encryptedID = crypt $encryptedID, $salt;
}
print OUT "sessionID=$encryptedID\n";
close OUT;
my $nLines1 = get_line_counts("$PATH_INFO/$loginname-config1.txt");
if($nLines1 && $nLines1 > $nLines*0.9){
	`cp $PATH_INFO/$loginname-config1.txt $PATH_INFO/$loginname-config.txt`;
}else{
	error_message("Failed to update configuration file!");
}
system("rm $PATH_INFO/$loginname-config1.txt");
return $sessionID;
}

#***********************************
sub  clean_up
#***********************************
{
my $date = `date \'+%y%m%d\'`;
chop($date);
my $date_cleaned;
if(open(INFO, "$PATH_INFO/cleaned_date.txt")){
	$date_cleaned=<INFO>;
	chomp $date_cleaned;
	close INFO;
}
if($date_cleaned==$date){ return; }
my @ls_info = `ls -l $PATH_INFO/guest`;
foreach my $filename (@ls_info){
	$filename =~ s/^.+guest/guest/;
	if($filename =~ /guest-config\.txt/){ next; }
	my $date_user = substr($filename,5,6);
	if($filename =~ /^guest/ && $date_user && $date > $date_user+1){
		`rm $PATH_INFO/$filename`;
	}
}
my @ls_data = `ls -l $PATH_DATA/guest`;
foreach my $filename (@ls_data){
	$filename =~ s/^.+guest/guest/;
	if($filename =~ /guest-config\.txt/){ next; }
	my $date_user = substr($filename,5,6);
	if($filename =~ /^guest/ && $date_user && $date > $date_user+1){
		`rm $PATH_DATA/$filename`;
	}
}
my $command = "ls -l --time-style=+\"%y%m%d\" $PATH_OUTPUT";
my @ls_output = `$command`;
foreach my $line (@ls_output){
	my @items = split(/\s+/,$line);
	my $n = @items;
	my $filename = $items[$n-1];
	my $filedate = $items[$n-2];
	if($filename =~ /^\d\d\d\d/ && $filedate < $date-1){
		`rm $PATH_OUTPUT/$filename`;
	}
}
`echo "$date" > $PATH_INFO/cleaned_date.txt`;
return;
}

#**************************************
sub  terminal_window
#**************************************
{
my $message = shift;
my $comment = shift;

print "<HTML><HEAD><TITLE>CisFinder response</TITLE>\n";
print_header();
print $message."\n";
print "<p><HR NOSHADE><p>\n";
print "<i>Please report any problems to <a href=mailto:sharov\@comcast.net>webmaster</a><p>\n";
if($comment eq "continue"){
	print "<FORM ACTION=$CGI_ADDRESS/cisfinder.cgi METHOD=POST>\n";
	print "<INPUT NAME=\"loginname\" TYPE=hidden VALUE=\"$loginname\">\n";
	print "<INPUT NAME=\"sessionID\" TYPE=hidden VALUE=\"$sessionID\">\n";
	print "<INPUT NAME=\"action\" TYPE=hidden VALUE=\"continue\">\n";
	print "<INPUT NAME=\"continue\" TYPE=submit VALUE=\"   Continue   \">\n";
	print "</FORM><p>\n";
}elsif($comment eq "register"){
	print "&nbsp; &nbsp; <A HREF=../index.html><IMG SRC=../images/button_home.gif BORDER=0></A><p>";
}else{
	print "<INPUT NAME=\"close_button\" TYPE=button VALUE=\"Close window\" style=width:160px; LANGUAGE=\"javascript\" onClick=\"window.close();\"><p>\n";
}
print "</BODY>\n";
print "</HTML>\n";
exit(0);
}

#**************************************
sub check_sessionID
#**************************************
{
my $loginname = shift;
my $sessionID = shift;

open(INFO, "$PATH_INFO/$loginname-config.txt") or error_message("Configuration file not found!");
my $savedID;
while(my $line = <INFO>){
	if($line =~ /^sessionID=/){
		$savedID = $line;
		$savedID =~ s/^sessionID=//;
		chomp $savedID;
		last;
	}
}
close INFO;
if(!$savedID){ error_message("SessionID not found!"); }
my $encryptedID = $sessionID;
for(my $i=0; $i<47; $i++){
	$encryptedID = crypt $encryptedID, $savedID;
}
if($encryptedID ne $savedID){
	error_message("Wrong sessionID! Only one session is allowed per user");
}
return(0);
}

#**********************************************************
sub  get_line_counts
#**********************************************************
{
my $filename = shift;
my $response = `wc -l $filename`;
$response =~ s/\D.+$//;
return $response;
}

