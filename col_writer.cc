#include "col_writer.hh"

/**
 * Constructs a col_writer instance of the templated number
 * of columns with the provided amount of padding
 * 
 * @param[in] padding the amount of padding for each column
 */
template <int N> col_writer<N>::col_writer(int padding_):dashes(NULL) {

	// initialize info for each column
	for (int i=0;i<N;i++) {
		padding[i]=padding_;
		headers[i]="";
		header_len[i]=0;
		max[i]=0;
	}

	setup();
}

/**
 * Class destructor
 */
template<int N> col_writer<N>::~col_writer() {
	if (dashes!=NULL) delete[] dashes;
}

/**
 * Sets up internal variables to prepare the col_writer
 * for outputting information. Can be run repeatedly
 * after updating headers or column max entry lengths
 */
template<int N> void col_writer<N>::setup() {

	// remove dash array before reallocating
	if (dashes!=NULL) delete[] dashes;

	// proceed through columns, calculating various
	// padding measures relative to headers/dashes
	for(int i=0;i<N;i++) {

		// length of dash underline should be the larger
		// of the header length and max entry length
		if (header_len[i]<=max[i]) {

			// in this case, need no padding inside dashes
			dash_len[i] = max[i];
			lr_dpad[i][0]=lr_dpad[i][1]=0;
		} else {

			// in this case, split "spare" spaces
			// between left and right of max entry
			dash_len[i] = header_len[i];
			lr_dpad[i][0] = (dash_len[i]-max[i])/2;
			lr_dpad[i][1] = dash_len[i]-max[i]-lr_dpad[i][0];
		}

		// split "spare" spaces in header line around left and right
		lr_hpad[i][0] = (max[i]-header_len[i])/2;
		lr_hpad[i][1] = max[i]-header_len[i]-lr_hpad[i][0];
	}

	// adjust padding until there's 2 spaces of padding on either
	// side of dash underline (i.e. 4 spaces between underlines)
	for (int i=0;i<N;i++) while (dlpad(i)<2) padding[i]++;

	// find max dash underline length
	int max_dashes=0;
	for(int i=0;i<N;i++) if(max_dashes<dash_len[i]) max_dashes=dash_len[i];

	// make C-string of that length
	// containing nothing but dashes
	dashes=new char[max_dashes+1];
	for(int i=0;i<max_dashes;i++) dashes[i]='-';
	dashes[max_dashes]='\0';
}

/**
 * Prints the set of N column headers
 *
 * @param[in] fh FILE stream to which headers are printed
 */
template<int N> void col_writer<N>::print_headers(FILE *fh) {

	// print each header
	for(int i=0;i<N;i++) {
		fprintf(fh,"%*s",hlpad(i),"");
		fprintf(fh,"%s",headers[i]);
	}
	fprintf(fh,"\n");

	// print each dash underline
	for(int i=0;i<N;i++) {
		fprintf(fh,"%*s",dlpad(i),"");
		fprintf(fh,"%.*s",dash_len[i],dashes);
	}
	fprintf(fh,"\n");
}

/**
 * Prints the set of N values in a given row
 *
 * @param[in] entries array of N C-strings representing the row
 * @param[in] fh FILE stream to which strings are printed
 */
template<int N> void col_writer<N>::print_row(
	char* entries[N],FILE *fh) {

	// print each string, adding padding
	// as necessary, and drop a new line
	for(int i=0;i<N;i++) {
		fprintf(fh,"%*s",padding[i],"");
		fprintf(fh,"%s%*s",entries[i],diff(entries[i],i),"");
	}
	fprintf(fh,"\n");
}

// explicit instantiation
template class col_writer<2>;
template class col_writer<3>;
template class col_writer<4>;
template class col_writer<5>;
template class col_writer<6>;
