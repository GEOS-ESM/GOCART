module FieldDictMap_mod
    use FieldDictEntry_mod
#include "types/key_deferredLengthString.inc"
#define _value class(FieldDictEntry)
#define _value_allocatable
#define _map FieldDict
#define _iterator FieldDictIterator
#define _alt
#include "templates/map.inc"
end module FieldDictMap_mod