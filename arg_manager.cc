///////////////////////////////////////////////////////////////////////////////
// arg_manager.cc
//
// The arg_manager class is used to read in and parse C++ command-line args
//
// Nick Derr, Nov 12, 2022
///////////////////////////////////////////////////////////////////////////////

#include "arg_manager.hh"
#include <ctype.h>

/**
 * Construct an instance of an arg_manager, given a count and
 * and values array of command-line args
 *
 * @param[in] c count of args
 * @param[in] v array of args
 */
arg_manager::arg_manager(int c,char** v,int n_req_):
	argc(c),n_req(n_req_),n_read(0),argv(v),used(new bool[argc]),
	help_ready(false),args_good(true),cw(2) {

	// initialize usage mask array
	for(int i=0;i<argc;i++)used[i]=false;

	// add help switches
	add_switch('h',"print help summary",
			"print a short listing of all args, "
			"types, and default values");
	add_switch('H',"print detailed help",
			"print a listing of all args, default values, "
			"and types, plus detailed descriptions");
	help_ready=true;

	// add column writer headers
	cw.set_header("flag",0);
	cw.set_header("name",1);
	cw.set_header("default",2);
	cw.set_header("type",3);

	// set length for flag column
	cw.update_col_len("-c",0);
}

/**
 * Class destructor
 */
arg_manager::~arg_manager() {

	// remove each arg allocation (for complex types)
	for (unsigned i=0;i<allocs.size();i++) delete[] allocs[i];

	// remove usage mask array
	delete[] used;
}

/**
 * Add a string as a required argument
 *
 * @param[in] name name of the argument
 */
void arg_manager::add_string(const char* name,const char* desc) {

	// check arg number
	if (n_read==n_req) {
		quit("error: attempt to add required arg #%d, but %d required args are"
				" specified\n",n_read+1,n_req);
	}

	// can't cast 0 to string so need to use string-specific add func
	char temp[MAX_STR_LEN];
	temp[0] = '\0';
	add_option<char*>(to_flag(n_read++),temp,name,desc);

	// this increments n_read!

	// generate a new C-string in memory and check the
	// submission doesn't overflow buffer
	allocs.push_back(new char[MAX_STR_LEN]);
	check_str(strings.back().val,to_flag(n_read-1)); // note n_read-1!

	// if all good, copy submitted string into new C-string
	// storage and point option at that
	sprintf(allocs.back(),"%s",strings.back().val);
	strings.back().val = allocs.back();
}

/**
 * Add a string arg with the provided flag, default value,
 * name, and description. This implementation accepts a
 * string literal as the default value
 *
 * @param[in] flag command-line flag for this option
 * @param[in] def default value if option not provided
 * @param[in] name name of the command-line option
 * @param[in] desc long-form (i.e. mult-line) description of option
 */
void arg_manager::add_string(char flag,const char *def,
		const char* name,const char* desc) {

	// opts need letter flags
	check_letter(flag);

	// alteration to hack compiler warning - accept string
	// literal for default value, but copy into a malleable
	// C-string before passing to the char* def version
	// automatically generated from templating
	char def_noconst[MAX_STR_LEN];

	// double-check the submitted default doesn't overflow buffer
	check_def(def,flag);
	sprintf(def_noconst,"%s",def);

	// use templated string submission function
	add_string(flag,def_noconst,name,desc);
}

/**
 * Add a string arg with the provided flag, default value,
 * name, and description. This implementation only accepts
 * a char* (i.e. not a string literal) as the default
 *
 * @param[in] flag command-line flag for this option
 * @param[in] def default value if option not provided
 * @param[in] name name of the command-line option
 * @param[in] desc long-form (i.e. mult-line) description of option
 */
void arg_manager::add_string(char flag,char* def,
		const char* name,const char* desc) {

	// opts need a letter flag
	check_letter(flag);

	// add the option
	add_option<char*>(flag,def,name,desc);

	// generate a new C-string in memory and check the
	// submission doesn't overflow buffer
	allocs.push_back(new char[MAX_STR_LEN]);
	check_str(strings.back().val,flag);

	// if all good, copy submitted string into new C-string
	// storage and point option at that
	sprintf(allocs.back(),"%s",strings.back().val);
	strings.back().val = allocs.back();
}

/**
 * Prints (with an fprintf-style format code) the
 * provided error message (and usage) and quits
 *
 * @param[in] format code
 * @param[in] ... additional arguments for format code
 */
void arg_manager::quit(const char * format,...) {

	printf("\n");

	va_list args;
	va_start(args,format);
	vfprintf(stderr,format,args);
	va_end(args);

	printf("\n\nHelp:\n\t%s -h    or    %s -H\n\n",
			argv[0],argv[0]);

	exit(-1);
}

/**
 * Checks a provided string to see if it is the provided
 * character in flag form (i.e. a hyphen followed by the char)
 *
 * @param[in] fl the character to test for
 * @param[in] str the string to test
 * @return whether the string is exactly "-<char>"
 */
inline bool arg_manager::is_flag(char fl,char* str) {
	char flag[3];
	sprintf(flag,"-%c",fl);
	return strcmp(flag,str)==0;
}

/**
 * Returns the index wi/in the arg list of the provided
 * flag, in the form "-<char>"
 *
 * @param[in] fl the provided flag
 * @return its index in the arg list
 */
int arg_manager::flag_pos(char fl) {

	// if digit, this is a required arg
	if (is_digit(fl)) return to_int(fl);

	// otherwise, look through list

	// proceed through arg list, keeping track of
	// whether or not we've already found the flag
	bool found=false;
	int pos=0;

	// (look for help flags everywhere, but otherwise only
	// look for flags in non-required slots)
	int i_max = (fl=='h'||fl=='H')?argc:(argc-n_req);
	for(int i=1; i<i_max; i++) {

		// if we find it, mark it! if we find it again, quit
		if(is_flag(fl,argv[i])) {
			if (!found) {
				pos=i;
				found=true;
			} else {
				quit("error: flag -%c provided multiple times\n",fl);
			}
		}
	}

	return pos;
}

/**
 * Returns the value (provided or default as appropriate)
 * associated with the option specifed by the provided flag
 *
 * @tparam C the type of the desired option
 * @param[in] fl the flag corresponding to the desired option
 * @return the value associated with the desired option
 */
template<class C> C arg_manager::option_value(char fl) {

	// grab pointer to option
	option_base *o = opt(fl);

	// make sure option type is consistent with function templating
	if (!o->is<C>()) {

		// error message depends on whether required
		if (!is_digit(fl)) {
			quit("error: call for type %s from flag -%c "
				"(option \"%s\"), which is defined as type %s\n",
				type<C>::name,fl,o->name,o->type_name());
		} else {
			quit("error: call for type %s from required arg #%d, "
				"which is defined as type %s\n",
				type<C>::name,to_int(fl)+1,o->type_name());
		}
	}

	// cast and return
	option<C> *o_cast = static_cast<option<C>*>(o);
	return o_cast->val;
}

/**
 * Returns a pointer to the option with the provided flag
 *
 * @param[in] flag the flag of the desired option
 * @return a pointer to the option associated with flag
 */
option_base* arg_manager::opt(char flag) {

	// use flag to look up type/index pair of option
	std::map<char,std::pair<type_ID,int> >::iterator ret;
	ret = opts.find(flag);

	// check that we found something
	if (ret==opts.end()) quit("error: no option with flag -%c found\n",flag);

	// return reference to option by indexing into vectors
	type_ID typ = ret->second.first;
	int ind  = ret->second.second;
	return opt(typ,ind);
}

/**
 * Returns a pointer to the option of the specified type with
 * the specified index in (type-specific) storage
 *
 * @param[in] type the type of this option
 * @param[in] ind the desired option's index in the specified
 *   type's storage
 * @return a pointer to the option at the specified index in
 *   the specified type's storage
 */
option_base* arg_manager::opt(type_ID typ,int ind) {

	// spit out location in appropriate vector
	switch (typ) {
		case DOUBLE: return &doubles[ind]; break;
		case INT:    return &ints[ind];    break;
		case STRING: return &strings[ind]; break;
		case BOOL: return &bools[ind]; break;
	}

	// print an error and quit if we get here
	quit("error: bad type enum");
	return NULL;
}

/**
 * Prints a message with usage information
 */
void arg_manager::print_usage() {

	// label and indent
	printf("\nUsage:\n\t");

	// function name and options
	printf("%s [options]",argv[0]);

	// required args
	for (int i=0;i<n_req;i++) {
		printf(" %s",opt(to_flag(i))->name);
	}

	// new-line
	printf("\n\n");

	// print arg info if there are any and they're all in
	if (n_req < 1 && n_read == n_req) return;


	printf("Args:\n");
	char str[MAX_STR_LEN];
	for (int i=0;i<n_req;i++) {

		// grab pointer to struct
		option_base *o = opt(to_flag(i));
		sprintf(str,"\t%s (%s): %s",o->name,o->type_name(),o->desc);
		print_indent_wrap(str,6,COLS);
	}
	printf("\n");
}

/**
 * Prints the help information for the arg_manager using the
 * provided names and descriptions. If the verbose flag is
 * provided, prints a long-form explanation for each option
 *
 * @param[in] verbose whether to print long-form descriptions
 */
void arg_manager::print_help(bool verbose) {

	print_usage();

	printf("Options:\n\n");

	// only print headers if we don't have
	// any long-form descriptions
	if(!verbose) cw.print_headers();

	// allocate strings for printing info
	char row_[4][MAX_STR_LEN];
	char* row[4] = {row_[0],row_[1],row_[2],row_[3]};
	char help[MAX_STR_LEN];

	// go through all options
	int fl;
	option_base *o;
	std::map<char,std::pair<type_ID,int> >::iterator it;
	for(it=opts.begin();it!=opts.end();it++) if (is_letter(it->first)) {

		// grab pointer to option
		fl = it->first;
		o = opt(fl);

		// store info to print in strings
		sprintf(row[0],"-%c",fl);
		sprintf(row[1],"%s",o->name);
		o->sprint_def(help,0);
		sprintf(row[2],"(default: %s)",help);
		sprintf(row[3],"%s",o->type_name());

		// print out nice columns
		cw.print_row(row,stdout);

		// if verbose, print wrapped description
		if (verbose) {
			printf("\n");
			print_indent_wrap(o->desc,INDENT,COLS);
			printf("\n");
		}
	}
}

/**
 * Processes arg to look for the help options
 */
void arg_manager::process_help() {

	// if help flag is provided, just print info and exit
	if (present('H')||present('h')) {
		print_help(present('H'));
		exit(0);
	}
}

/**
 * Processes all submitted arguments, making sure that each
 * corresponds to a valid option. Takes any "short-circuit"
 * actions (e.g. printing help info) if requested. If this
 * runs and the program does not display an error, it can
 * be assumed that any provided arguments were valid and that
 * the retrievable answer of each option corresponds to
 * either a user-submitted value or its default
 */
void arg_manager::process_args() {

	// check we got the right number of required args
	if (n_read != n_req) {
		quit("error: %d required args specified, but only %d provided\n",
				n_req,n_read);
	}

	// set up column writer
	cw.setup();

	// check for help
	process_help();

	// otherwise, check to see that every submitted command-line
	// argument corresponds to a valid option
	bool all_used = true;
	for(int i=1;i<argc-n_req;i++) if (!used[i]) {
		all_used = false;
		if (i<argc-n_req) {
			if (strlen(argv[i])==2 && argv[i][0]=='-') {
				printf("unrecognized flag:     %s\n",argv[i]);
			} else {
				printf("provided without flag: %s\n",argv[i]);
			}
		}
	}

	// bail if any weird args
	if (!all_used) quit("error: unused or unrecognized args provided\n");

	// now check required args
	process_reqs();
}

/**
 * Processes required (positional) arguments.
 */
void arg_manager::process_reqs() {

	if (!args_good) {
		print_usage();
		quit("error: incorrect number of provided args\n");
	}

	// array for type name
	char tname[MAX_STR_LEN];

	for (int i=0;i<n_req;i++) {

		// get pointer to struct and index in arg array
		option_base *o = opt(to_flag(i));
		int ai = argc-n_req+i;

		bool success;

		if (o->is<double>()) {
			option<double> *opt=static_cast<option<double>*>(o);
			success = opt->parse(argv[ai]);
			sprintf(tname,"%s",opt->type_name());
		} else if (o->is<int>()) {
			option<int> *opt=static_cast<option<int>*>(o);
			success = opt->parse(argv[ai]);
			sprintf(tname,"%s",opt->type_name());
		} else if (o->is<char*>()) {
			option<char*> *opt=static_cast<option<char*>*>(o);
			success = opt->parse(argv[ai]);
			sprintf(tname,"%s",opt->type_name());
		} else {
			quit("error: bad required arg type\n");
		}

		if (!success) {
			quit("error: required arg #%d (\"%s\") has provided value"
				" \"%s\", which cannot be parsed as a %s\n",
				i+1,o->name,argv[ai],tname);
		}
	}
}

/**
 * Prints the provided C-string, indented the provided number
 * of spaces and wrapped (to the nearest whitespace) after
 * the provided number of columns. Any newlines within the string
 * are printed as newlines within the output.
 *
 * @param[in] str the C-string to print
 * @param[in] indent the number of columns to indent the output
 * @param[in] cols the number of columns after which to wrap
 */
void arg_manager::print_indent_wrap(const char *str,unsigned indent,unsigned cols) {

	// grab length of output, variable marking length
	// until newline, and pointer with which to proceed
	// through string
	size_t len = strlen(str);
	const char* ss = str;
	int line = 0;

	do {

		// move pointer until at first non-
		// whitespace or newline character
		while (isspace(*ss) && *ss!='\0'
				&& !(ss!=str&&ss[-1]=='\n')) ss++;

		// get distance to end of string
		size_t to_go = strlen(ss);

		// if we're not at the end, figure
		// out how long the string should be
		if (to_go!=0) {

			if (to_go < cols) {

				// if we can avoid wrapping, we're done
				line = to_go;
			} else {

				// otherwise, pick the first non-whitespace
				// bit before the wrapping limit and wrap it
				line = cols;
				while (line>0 && ss[line]!=' ') line--;
				if (line==0) line = cols;
			}

			// look for any included newlines
			for (int i=0;i<line;i++) if (ss[i]=='\n') {
				line=i;
				break;
			}

			// print the selected portion with a newline
			printf("%*c%.*s\n",indent,' ',line,ss);
			ss += line;
		}

	} while (ss<str+len);
}

/**
 * Returns a pointer to the vector containing options of the
 * provided type
 *
 * @tparam T the provided type
 * @return a pointer to the vector containing T options
 */
template<class C> std::vector<option<C> >* arg_manager::opt_storage() {

	// type aliases for storage vectors
	typedef std::vector<option<C> > vec;
	typedef std::vector<option<C> > vec;

	// returns casts based on type ID (this will
	// resolve at compile time)
	switch (type<C>::ID) {
		case INT   : return static_cast<vec*>(static_cast<void*>(&ints));
		case DOUBLE: return static_cast<vec*>(static_cast<void*>(&doubles));
		case STRING: return static_cast<vec*>(static_cast<void*>(&strings));
		case BOOL  : return static_cast<vec*>(static_cast<void*>(&bools));
	}

	// squash compilation warning with spuriours return
	quit("error: bad type template parameter\n");
	return NULL;
}

/**
 * Checks a provided string to see if its length is less than the
 * maxmimum allowed with the current C-string implementation. The
 * error message which prints specifies any violators as provided
 * or default according to the provided argument
 *
 * @param[in] str the string to check
 * @param[in] flag the flag associated with the string in question
 * @param[in] def whether the value being tested is a default (default=false)
 */
inline void arg_manager::check_str(const char *str,char flag,bool def) {
	if (strlen(str) < MAX_STR_LEN) return;
	quit("error: %s value for string arg -%c exceeds max length %d\n",
		def?"default":"provided",flag,MAX_STR_LEN);
}

/**
 * Bool-specific implementation of option_set, which marks an option as
 * having been submitted by the user and its corresponding arguments
 * as having been used
 *
 * @param[in] a the option to mark the submission of
 */
template<> void arg_manager::option_set<bool>(option<bool> &a) {

	// grab option's place in command-line args
	int pos = flag_pos(a.flag);

	// if this is 0, it wasn't in the args
	// so we don't need to mark as used
	if (pos==0) return;

	// if it's here, mark it
	a.present = true;

	// a switch becomes true if it's here
	a.val = true;

	// mark the arg as having been used
	used[pos] = true;
}

/**
 * Marks an option as having been submitted by the user, parses
 * the submitted string and assigns its value to the relevent option,
 * and marks the corresponding arguments as having been used
 *
 * @param[in] a the option to mark the submission of
 * @param[in] req whether this is a required argument
 */
template<class C> void arg_manager::option_set(option<C> &a) {

	// grab option's place in command-line args and verify
	// that it has actually been submitted
	int pos = flag_pos(a.flag);

	// if this is 0, it wasn't in the args
	// so we don't need to mark as used
	if (pos==0) return;

	// if it's here, mark it
	a.present = true;

	// check that it's got an associated value with it
	if(++pos==argc-n_req) quit("error: providing the flag -%c for arg \"%s\" "
				"of type %s requires an additional argument\n",
				a.flag,a.name,typeid(C).name());

	// if help not present, check that the value can be parsed
	if(!a.parse(argv[pos])) {

		// error mesg depends if required
		quit("error: in reading in arg \"%s\""
			" with flag -%c, failed to parse string %s as"
			" %s\n",a.name,a.flag,argv[pos],a.type_name());
	}

	// mark flag and additional arg as used
	used[pos-1]=true;
	used[pos]=true;
}

/**
 * Adds an option of the specified type to the arg_manager, such that
 * it has the provided flag, default value, name, and long-form
 * description
 *
 * @param[in] flag the option's flag (must be unique)
 * @param[in] def the option's default value
 * @param[in] name the option's name
 * @param[in] descr an long-form description of the option's function
 */
template<class C> void arg_manager::add_option(
		char flag,C def,const char* name,const char* descr) {

	// grab a pointer to the storage vector for this datatype
	// and try to insert an option with this flag
	std::vector<option<C> > *v = opt_storage<C>();
	std::pair<std::map<char,std::pair<type_ID,int> >::iterator,bool> ret;
	ret = opts.insert(std::pair<char,std::pair<type_ID,int> >(flag,
				std::pair<type_ID,int>(type<C>::ID,v->size())));

	// check that the insertion was successful. if not, this
	// flag has already been used. inform user and quit
	if (!ret.second) {
		type_ID type = ret.first->second.first;
		int ind  = ret.first->second.second;

		// grab option already used with flag for info
		option_base *o = opt(type,ind);
		quit("error: flag -%c already used for option \"%s\" of type %s\n",
				flag,o->name,o->type_name());
	}

	// if flag is not a duplicate, store the option and mark its submission
	v->push_back(option<C>(flag,def,name,descr));
	option<C>&o = v->back();

	// if this is a required arg (flag is digit) check index in range
	// otherwise mark optional arg as having been submitted
	if(is_digit(flag)) {
		int pos=argc-n_req+to_int(flag);
		if(pos<1||pos>=argc) args_good=false;
	} else {
		option_set<C>(o);
	}

	// record lengths for nice printing
	char dval[MAX_STR_LEN],def_str[MAX_STR_LEN];
	o.sprint_def(dval,0);
	sprintf(def_str,"(default: %s)",dval);

	if (is_letter(flag)) {
		cw.update_col_len(o.name,1);
		cw.update_col_len(def_str,2);
		cw.update_col_len(o.type_name(),3);
	}
}

// explicit instantiations
template double arg_manager::option_value<double>(char fl);
template int arg_manager::option_value<int>(char fl);
template char* arg_manager::option_value<char*>(char fl);
template bool arg_manager::option_value<bool>(char fl);

template void arg_manager::add_option<double>(char,double,const char*,const char*);
template void arg_manager::add_option<int>(char,int,const char*,const char*);
template void arg_manager::add_option<char*>(char,char*,const char*,const char*);
template void arg_manager::add_option<bool>(char,bool,const char*,const char*);
