#include <locale.h>

#include "../src/seqevtpair.h"

/* completely nonsense numbers, just used to test I/O */
const char* xmlIn = "<?xml version=\"1.0\"?>\n"
"<model><order>2</order>"
"<delete><begin>0.1</begin><extend>0.9</extend></delete>"
"<states>"
"<state><id>AA</id><emit>0.11</emit><mean>1</mean><precision>2</precision></state>"
"<state><id>AC</id><emit>0.12</emit><mean>3</mean><precision>21</precision></state>"
"<state><id>AG</id><emit>0.13</emit><mean>5</mean><precision>22</precision></state>"
"<state><id>AT</id><emit>0.14</emit><mean>7</mean><precision>23</precision></state>"
"<state><id>CA</id><emit>0.15</emit><mean>9</mean><precision>24</precision></state>"
"<state><id>CC</id><emit>0.16</emit><mean>11</mean><precision>25</precision></state>"
"<state><id>CG</id><emit>0.17</emit><mean>13</mean><precision>26</precision></state>"
"<state><id>CT</id><emit>0.18</emit><mean>15</mean><precision>27</precision></state>"
"<state><id>GA</id><emit>0.19</emit><mean>17</mean><precision>28</precision></state>"
"<state><id>GC</id><emit>0.21</emit><mean>19</mean><precision>29</precision></state>"
"<state><id>GG</id><emit>0.22</emit><mean>21</mean><precision>32</precision></state>"
"<state><id>GT</id><emit>0.23</emit><mean>31</mean><precision>42</precision></state>"
"<state><id>TA</id><emit>0.24</emit><mean>41</mean><precision>52</precision></state>"
"<state><id>TC</id><emit>0.25</emit><mean>51</mean><precision>62</precision></state>"
"<state><id>TG</id><emit>0.26</emit><mean>61</mean><precision>72</precision></state>"
"<state><id>TT</id><emit>0.27</emit><mean>71</mean><precision>82</precision></state>"
"</states>"
"<start><emit>0.28</emit><mean>72</mean><precision>88</precision></start>"
"</model>\n";

int main(int argc, char * argv[])
{
  Seq_event_pair_model* model;
  xmlChar *xmlOut;
  model = new_seq_event_pair_model_from_xml_string (xmlIn);
  xmlOut = convert_seq_event_pair_model_to_xml_string (model);
  if (strcmp((char*)xmlIn,(char*)xmlOut) == 0) {
    printf ("ok\n");
  } else {
    printf ("not ok\n");
    fprintf (stderr, "Expected output:\n%s\n\nActual output:\n%s\n", xmlIn, xmlOut);
  }
  SafeFree (xmlOut);
  delete_seq_event_pair_model (model);
  return EXIT_SUCCESS;
}
