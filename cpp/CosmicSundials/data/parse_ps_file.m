%{

file: parse_ps_file.m
author: 2012 Ben O'Hara
contact: buohara@gmail.com, buohara@uvm.edu
description: Parses a ps file and writes its contents to a separate .txt file in a format more conducive for use in C.

@param ps_name Name of the power system to parse. This name is the same as the function name in whichever ps file is used (e.g., function ps = ***case39***), and as such, the ps file should exist in the same directory as this file.

@param ps_file_out Name of the file where output will be placed.

%}

function parse_ps_file(ps_name, ps_file_out)

	%try and load the ps struct
	ps_handle = str2func(ps_name);
	ps = feval(ps_handle);
	
	%grab the base MVA
	base_mva = ps.baseMVA;
	
	%grab the sizes of bus structs, branch structs, etc.
	num_buses = length(ps.bus(:,1));
	num_branches = length(ps.branch(:,1));
	num_gens = length(ps.gen(:,1));

	%not all structs exist in all files, so handle entries that may or may not exist.
	
	%mac
	macs_found = false;
	try
		num_macs = length(ps.mac(:,1));
		macs_found = true;
	catch exception
		fprintf('Warning: Unable to locate machine data in file %s.m.\n', ps_name);
	end
	
	%{
	This script operates under the assumption that shunt data is available either in its own substructure (ps.shunt) or is embedded in the bus struct. The following shunts_available flag indicates whether the data is available (true) or if we have to get it from buses (false).
	%}
	shunts_found = false;
	try
		num_shunts = length(ps.shunt(:,1));
		shunts_available = true;
	catch exception
		fprintf('Warning: Unable to locate shunt data in file %s.m. Using Pd and Qd data from ps.bus.\n', ps_name);
	end
	
	%if shunt data is unavailable, extract it from bus data
	if ~shunts_found
		shunts = ps.bus(ps.bus(:,3)~=0 | ps.bus(:,4)~=0,[1 3 4]);
		num_shunts = length(shunts(:,1));
	end

	%now begin writing output to file
	out_file = fopen(ps_file_out, 'w');
	
	%print base MVA
	fprintf(out_file, 'BASE_MVA %g', base_mva);
	
	%print counts of each type of entry
	fprintf(out_file, 'BUS %d\n', num_buses);
	fprintf(out_file, 'BRANCH %d\n', num_branches);
	fprintf(out_file, 'GEN %d\n', num_gens);
	fprintf(out_file, 'SHUNT %d\n', num_shunts);
	
	fprintf(out_file, '\n');
	
	%now print data to the file
	fprintf(out_file, '%d %d %4.5f %4.5f %4.5f %4.5f %d %4.5f %4.5f %4.5f %d %4.5f %4.5f\n', ps.bus(:,1:1:13)');
	fprintf(out_file, '\n');
	
	fprintf(out_file, '%d %d %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f\n',ps.branch(:,1:1:8)');
	fprintf(out_file, '\n');
	
	fprintf(out_file, '%d %4.5f %4.5f %4.5f %4.5f %4.5f %4.5f %d %4.5f %4.5f\n', ps.gen(:,1:1:10)');
	fprintf(out_file, '\n');
	
	if shunts_found
		fprintf(out_file, '%d %4.5f %4.5f\n', ps.shunt(:,1:1:3)');
	else
		fprintf(out_file, '%d %4.5f %4.5f\n', shunts(:,1:1:3)');
	end
end