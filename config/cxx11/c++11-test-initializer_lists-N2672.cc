#include <map>
#include <string>
#include <vector>
#include <utility>
#include <iostream>

using namespace std;

int main()
{
    bool allok(true);
    
    map<string, vector<pair<string, int>>> name_languages_year {
        {"Dennis Ritchie",    {{"B",      1969}, {"C",        1973}}},
        {"Niklaus Wirth",     {{"Pascal", 1970}, {"Modula-2", 1973}, {"Oberon", 1986}}},
        {"Bjarne Stroustrup", {{"C++",    1983}}},
        {"Walter Bright",     {{"D",      1999}}}
    };
    // notice how the lists are nested to match the templates' parameters

    // adds a new entry to the map
    name_languages_year["John McCarthy"] = {
        {"Lisp", 1958}
    };
    // notice the lack of explicit types

    if( name_languages_year["Niklaus Wirth"].at(0).first != "Pascal" )
        allok=false;
    if( name_languages_year["Niklaus Wirth"].at(1).first != "Modula-2" )
        allok=false;
    if( name_languages_year["Niklaus Wirth"].at(2).first != "Oberon" )
        allok=false;
    if( name_languages_year["John McCarthy"].at(0).first != "Lisp" )
        allok=false;
    if( name_languages_year["John McCarthy"].at(0).second != 1958 )
        allok=false;

    {
        std::vector<std::string> v1 = { "xyzzy", "plugh", "abracadabra" };
        std::vector<std::string> v2({ "xyzzy", "plugh", "abracadabra" });
        std::vector<std::string> v3{ "xyzzy", "plugh", "abracadabra" };

        if( v1 != v2 || v1 != v3 || v2 != v3) allok=false;

        // range-based for
        for( std::string &item: v1 )
            std::cout << item << endl;
    }
    
    return allok ? 0 : 1;
}
