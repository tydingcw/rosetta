Ř
Î##
:
Add
x"T
y"T
z"T"
Ttype:
2	

ArgMax

input"T
	dimension"Tidx
output"output_type" 
Ttype:
2	"
Tidxtype0:
2	"
output_typetype0	:
2	
P
Assert
	condition
	
data2T"
T
list(type)(0"
	summarizeint
E
AssignAddVariableOp
resource
value"dtype"
dtypetype
B
AssignVariableOp
resource
value"dtype"
dtypetype
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
~
BiasAddGrad
out_backprop"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
R
BroadcastGradientArgs
s0"T
s1"T
r0"T
r1"T"
Ttype0:
2	
N
Cast	
x"SrcT	
y"DstT"
SrcTtype"
DstTtype"
Truncatebool( 
h
ConcatV2
values"T*N
axis"Tidx
output"T"
Nint(0"	
Ttype"
Tidxtype0:
2	
8
Const
output"dtype"
valuetensor"
dtypetype
š
DenseToDenseSetOperation	
set1"T	
set2"T
result_indices	
result_values"T
result_shape	"
set_operationstring"
validate_indicesbool("
Ttype:
	2	
5
DivNoNan
x"T
y"T
z"T"
Ttype:
2
S
DynamicStitch
indices*N
data"T*N
merged"T"
Nint(0"	
Ttype
B
Equal
x"T
y"T
z
"
Ttype:
2	

W

ExpandDims

input"T
dim"Tdim
output"T"	
Ttype"
Tdimtype0:
2	
^
Fill
dims"
index_type

value"T
output"T"	
Ttype"

index_typetype0:
2	
?
FloorDiv
x"T
y"T
z"T"
Ttype:
2	
9
FloorMod
x"T
y"T
z"T"
Ttype:

2	
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
;
Maximum
x"T
y"T
z"T"
Ttype:

2	

Mean

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
N
Merge
inputs"T*N
output"T
value_index"	
Ttype"
Nint(0
=
Mul
x"T
y"T
z"T"
Ttype:
2	
.
Neg
x"T
y"T"
Ttype:

2	

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
X
PlaceholderWithDefault
input"dtype
output"dtype"
dtypetype"
shapeshape
6
Pow
x"T
y"T
z"T"
Ttype:

2	

Prod

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
~
RandomUniform

shape"T
output"dtype"
seedint "
seed2int "
dtypetype:
2"
Ttype:
2	
a
Range
start"Tidx
limit"Tidx
delta"Tidx
output"Tidx"
Tidxtype0:	
2	
@
ReadVariableOp
resource
value"dtype"
dtypetype
>
RealDiv
x"T
y"T
z"T"
Ttype:
2	
E
Relu
features"T
activations"T"
Ttype:
2	
V
ReluGrad
	gradients"T
features"T
	backprops"T"
Ttype:
2	
[
Reshape
tensor"T
shape"Tshape
output"T"	
Ttype"
Tshapetype0:
2	
ŕ
ResourceApplyAdam
var
m
v
beta1_power"T
beta2_power"T
lr"T

beta1"T

beta2"T
epsilon"T	
grad"T" 
Ttype:
2	"
use_lockingbool( "
use_nesterovbool( 
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
P
Shape

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
O
Size

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
G
SquaredDifference
x"T
y"T
z"T"
Ttype:

2	
:
Sub
x"T
y"T
z"T"
Ttype:
2	

Sum

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
M
Switch	
data"T
pred

output_false"T
output_true"T"	
Ttype
c
Tile

input"T
	multiples"
Tmultiples
output"T"	
Ttype"

Tmultiplestype0:
2	
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape
9
VarIsInitializedOp
resource
is_initialized
"train*2.0.0-alpha02v1.12.0-9492-g2c319fb4158Ř­
n
dense_inputPlaceholder*
dtype0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
*
shape:˙˙˙˙˙˙˙˙˙


-dense/kernel/Initializer/random_uniform/shapeConst*
valueB"
      *
_class
loc:@dense/kernel*
dtype0*
_output_shapes
:

+dense/kernel/Initializer/random_uniform/minConst*
valueB
 *č!ż*
_class
loc:@dense/kernel*
dtype0*
_output_shapes
: 

+dense/kernel/Initializer/random_uniform/maxConst*
_output_shapes
: *
valueB
 *č!?*
_class
loc:@dense/kernel*
dtype0
Ě
5dense/kernel/Initializer/random_uniform/RandomUniformRandomUniform-dense/kernel/Initializer/random_uniform/shape*
T0*
_class
loc:@dense/kernel*
dtype0*
_output_shapes

:

Î
+dense/kernel/Initializer/random_uniform/subSub+dense/kernel/Initializer/random_uniform/max+dense/kernel/Initializer/random_uniform/min*
_output_shapes
: *
T0*
_class
loc:@dense/kernel
ŕ
+dense/kernel/Initializer/random_uniform/mulMul5dense/kernel/Initializer/random_uniform/RandomUniform+dense/kernel/Initializer/random_uniform/sub*
T0*
_class
loc:@dense/kernel*
_output_shapes

:

Ň
'dense/kernel/Initializer/random_uniformAdd+dense/kernel/Initializer/random_uniform/mul+dense/kernel/Initializer/random_uniform/min*
T0*
_class
loc:@dense/kernel*
_output_shapes

:


dense/kernelVarHandleOp*
shape
:
*
shared_namedense/kernel*
_class
loc:@dense/kernel*
dtype0*
_output_shapes
: 
i
-dense/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense/kernel*
_output_shapes
: 

dense/kernel/AssignAssignVariableOpdense/kernel'dense/kernel/Initializer/random_uniform*
_class
loc:@dense/kernel*
dtype0

 dense/kernel/Read/ReadVariableOpReadVariableOpdense/kernel*
_output_shapes

:
*
_class
loc:@dense/kernel*
dtype0

dense/bias/Initializer/zerosConst*
valueB*    *
_class
loc:@dense/bias*
dtype0*
_output_shapes
:


dense/biasVarHandleOp*
shape:*
shared_name
dense/bias*
_class
loc:@dense/bias*
dtype0*
_output_shapes
: 
e
+dense/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOp
dense/bias*
_output_shapes
: 
{
dense/bias/AssignAssignVariableOp
dense/biasdense/bias/Initializer/zeros*
dtype0*
_class
loc:@dense/bias

dense/bias/Read/ReadVariableOpReadVariableOp
dense/bias*
_class
loc:@dense/bias*
dtype0*
_output_shapes
:
h
dense/MatMul/ReadVariableOpReadVariableOpdense/kernel*
dtype0*
_output_shapes

:

r
dense/MatMulMatMuldense_inputdense/MatMul/ReadVariableOp*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
c
dense/BiasAdd/ReadVariableOpReadVariableOp
dense/bias*
dtype0*
_output_shapes
:
v
dense/BiasAddBiasAdddense/MatMuldense/BiasAdd/ReadVariableOp*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
S

dense/ReluReludense/BiasAdd*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
Ł
/dense_1/kernel/Initializer/random_uniform/shapeConst*
valueB"      *!
_class
loc:@dense_1/kernel*
dtype0*
_output_shapes
:

-dense_1/kernel/Initializer/random_uniform/minConst*
valueB
 *mż*!
_class
loc:@dense_1/kernel*
dtype0*
_output_shapes
: 

-dense_1/kernel/Initializer/random_uniform/maxConst*
_output_shapes
: *
valueB
 *m?*!
_class
loc:@dense_1/kernel*
dtype0
Ň
7dense_1/kernel/Initializer/random_uniform/RandomUniformRandomUniform/dense_1/kernel/Initializer/random_uniform/shape*
T0*!
_class
loc:@dense_1/kernel*
dtype0*
_output_shapes

:
Ö
-dense_1/kernel/Initializer/random_uniform/subSub-dense_1/kernel/Initializer/random_uniform/max-dense_1/kernel/Initializer/random_uniform/min*
T0*!
_class
loc:@dense_1/kernel*
_output_shapes
: 
č
-dense_1/kernel/Initializer/random_uniform/mulMul7dense_1/kernel/Initializer/random_uniform/RandomUniform-dense_1/kernel/Initializer/random_uniform/sub*
T0*!
_class
loc:@dense_1/kernel*
_output_shapes

:
Ú
)dense_1/kernel/Initializer/random_uniformAdd-dense_1/kernel/Initializer/random_uniform/mul-dense_1/kernel/Initializer/random_uniform/min*
_output_shapes

:*
T0*!
_class
loc:@dense_1/kernel

dense_1/kernelVarHandleOp*
_output_shapes
: *
shape
:*
shared_namedense_1/kernel*!
_class
loc:@dense_1/kernel*
dtype0
m
/dense_1/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_1/kernel*
_output_shapes
: 

dense_1/kernel/AssignAssignVariableOpdense_1/kernel)dense_1/kernel/Initializer/random_uniform*!
_class
loc:@dense_1/kernel*
dtype0

"dense_1/kernel/Read/ReadVariableOpReadVariableOpdense_1/kernel*!
_class
loc:@dense_1/kernel*
dtype0*
_output_shapes

:

dense_1/bias/Initializer/zerosConst*
valueB*    *
_class
loc:@dense_1/bias*
dtype0*
_output_shapes
:

dense_1/biasVarHandleOp*
_output_shapes
: *
shape:*
shared_namedense_1/bias*
_class
loc:@dense_1/bias*
dtype0
i
-dense_1/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_1/bias*
_output_shapes
: 

dense_1/bias/AssignAssignVariableOpdense_1/biasdense_1/bias/Initializer/zeros*
_class
loc:@dense_1/bias*
dtype0

 dense_1/bias/Read/ReadVariableOpReadVariableOpdense_1/bias*
dtype0*
_output_shapes
:*
_class
loc:@dense_1/bias
l
dense_1/MatMul/ReadVariableOpReadVariableOpdense_1/kernel*
dtype0*
_output_shapes

:
u
dense_1/MatMulMatMul
dense/Reludense_1/MatMul/ReadVariableOp*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
g
dense_1/BiasAdd/ReadVariableOpReadVariableOpdense_1/bias*
dtype0*
_output_shapes
:
|
dense_1/BiasAddBiasAdddense_1/MatMuldense_1/BiasAdd/ReadVariableOp*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0

dense_1_targetPlaceholder*0
_output_shapes
:˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙*%
shape:˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙*
dtype0
R
ConstConst*
valueB*  ?*
dtype0*
_output_shapes
:

dense_1_sample_weightsPlaceholderWithDefaultConst*
dtype0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*
shape:˙˙˙˙˙˙˙˙˙
v
total/Initializer/zerosConst*
valueB
 *    *
_class

loc:@total*
dtype0*
_output_shapes
: 
x
totalVarHandleOp*
shared_nametotal*
_class

loc:@total*
dtype0*
_output_shapes
: *
shape: 
[
&total/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal*
_output_shapes
: 
g
total/AssignAssignVariableOptotaltotal/Initializer/zeros*
_class

loc:@total*
dtype0
q
total/Read/ReadVariableOpReadVariableOptotal*
_class

loc:@total*
dtype0*
_output_shapes
: 
v
count/Initializer/zerosConst*
valueB
 *    *
_class

loc:@count*
dtype0*
_output_shapes
: 
x
countVarHandleOp*
shared_namecount*
_class

loc:@count*
dtype0*
_output_shapes
: *
shape: 
[
&count/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount*
_output_shapes
: 
g
count/AssignAssignVariableOpcountcount/Initializer/zeros*
dtype0*
_class

loc:@count
q
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
_class

loc:@count*
dtype0
l
!metrics/accuracy/ArgMax/dimensionConst*
_output_shapes
: *
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0

metrics/accuracy/ArgMaxArgMaxdense_1_target!metrics/accuracy/ArgMax/dimension*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
n
#metrics/accuracy/ArgMax_1/dimensionConst*
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0*
_output_shapes
: 

metrics/accuracy/ArgMax_1ArgMaxdense_1/BiasAdd#metrics/accuracy/ArgMax_1/dimension*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙

metrics/accuracy/EqualEqualmetrics/accuracy/ArgMaxmetrics/accuracy/ArgMax_1*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0	
r
metrics/accuracy/CastCastmetrics/accuracy/Equal*

SrcT0
*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*

DstT0
`
metrics/accuracy/ConstConst*
valueB: *
dtype0*
_output_shapes
:
k
metrics/accuracy/SumSummetrics/accuracy/Castmetrics/accuracy/Const*
T0*
_output_shapes
: 
e
$metrics/accuracy/AssignAddVariableOpAssignAddVariableOptotalmetrics/accuracy/Sum*
dtype0

metrics/accuracy/ReadVariableOpReadVariableOptotal%^metrics/accuracy/AssignAddVariableOp^metrics/accuracy/Sum*
dtype0*
_output_shapes
: 
U
metrics/accuracy/SizeSizemetrics/accuracy/Cast*
T0*
_output_shapes
: 
f
metrics/accuracy/Cast_1Castmetrics/accuracy/Size*
_output_shapes
: *

DstT0*

SrcT0

&metrics/accuracy/AssignAddVariableOp_1AssignAddVariableOpcountmetrics/accuracy/Cast_1%^metrics/accuracy/AssignAddVariableOp*
dtype0
Ż
!metrics/accuracy/ReadVariableOp_1ReadVariableOpcount%^metrics/accuracy/AssignAddVariableOp'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

*metrics/accuracy/div_no_nan/ReadVariableOpReadVariableOptotal'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

,metrics/accuracy/div_no_nan/ReadVariableOp_1ReadVariableOpcount'^metrics/accuracy/AssignAddVariableOp_1*
_output_shapes
: *
dtype0
˘
metrics/accuracy/div_no_nanDivNoNan*metrics/accuracy/div_no_nan/ReadVariableOp,metrics/accuracy/div_no_nan/ReadVariableOp_1*
T0*
_output_shapes
: 
c
metrics/accuracy/IdentityIdentitymetrics/accuracy/div_no_nan*
_output_shapes
: *
T0
n
#metrics/accuracy/ArgMax_2/dimensionConst*
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0*
_output_shapes
: 

metrics/accuracy/ArgMax_2ArgMaxdense_1_target#metrics/accuracy/ArgMax_2/dimension*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
n
#metrics/accuracy/ArgMax_3/dimensionConst*
_output_shapes
: *
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0

metrics/accuracy/ArgMax_3ArgMaxdense_1/BiasAdd#metrics/accuracy/ArgMax_3/dimension*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0

metrics/accuracy/Equal_1Equalmetrics/accuracy/ArgMax_2metrics/accuracy/ArgMax_3*
T0	*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
v
metrics/accuracy/Cast_2Castmetrics/accuracy/Equal_1*

SrcT0
*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*

DstT0
b
metrics/accuracy/Const_1Const*
valueB: *
dtype0*
_output_shapes
:
q
metrics/accuracy/MeanMeanmetrics/accuracy/Cast_2metrics/accuracy/Const_1*
T0*
_output_shapes
: 

6loss/dense_1_loss/mean_squared_error/SquaredDifferenceSquaredDifferencedense_1/BiasAdddense_1_target*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0

;loss/dense_1_loss/mean_squared_error/Mean/reduction_indicesConst*
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0*
_output_shapes
: 
Ô
)loss/dense_1_loss/mean_squared_error/MeanMean6loss/dense_1_loss/mean_squared_error/SquaredDifference;loss/dense_1_loss/mean_squared_error/Mean/reduction_indices*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
­
gloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapeShapedense_1_sample_weights*
T0*
_output_shapes
:
¨
floss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rankConst*
dtype0*
_output_shapes
: *
value	B :
ż
floss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapeShape)loss/dense_1_loss/mean_squared_error/Mean*
_output_shapes
:*
T0
§
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rankConst*
dtype0*
_output_shapes
: *
value	B :
§
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar/xConst*
value	B : *
dtype0*
_output_shapes
: 
Ü
closs/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarEqualeloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar/xfloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank*
T0*
_output_shapes
: 
ć
oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/SwitchSwitchcloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarcloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar*
T0
*
_output_shapes
: : 

qloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_tIdentityqloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch:1*
T0
*
_output_shapes
: 

qloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_fIdentityoloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch*
_output_shapes
: *
T0


ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_idIdentitycloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar*
_output_shapes
: *
T0

í
qloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1Switchcloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
_output_shapes
: : *
T0
*v
_classl
jhloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar
î
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankEqualloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switchloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1*
_output_shapes
: *
T0

loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/SwitchSwitcheloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rankploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*x
_classn
ljloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank*
_output_shapes
: : *
T0

loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1Switchfloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rankploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
T0*y
_classo
mkloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank*
_output_shapes
: : 
Ű
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/SwitchSwitchloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
T0
*
_output_shapes
: : 
Ç
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_tIdentityloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:1*
T0
*
_output_shapes
: 
Ĺ
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_fIdentityloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch*
T0
*
_output_shapes
: 
Ę
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_idIdentityloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
T0
*
_output_shapes
: 
ý
˘loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dimConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0*
_output_shapes
: 
§
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims
ExpandDimsŠloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1˘loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim*
_output_shapes

:*
T0
°
Ľloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/SwitchSwitchfloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapeploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
T0*y
_classo
mkloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::

§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1SwitchĽloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switchloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0*y
_classo
mkloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::

Łloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/ShapeConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
_output_shapes
:*
valueB"      *
dtype0
ő
Łloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/ConstConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
value	B :*
dtype0*
_output_shapes
: 

loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_likeFillŁloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/ShapeŁloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const*
T0*
_output_shapes

:
ń
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axisConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
value	B :*
dtype0*
_output_shapes
: 
ź
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concatConcatV2loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDimsloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_likeloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis*
T0*
N*
_output_shapes

:
˙
¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dimConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0*
_output_shapes
: 
­
 loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1
ExpandDimsŤloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim*
T0*
_output_shapes

:
´
§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/SwitchSwitchgloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapeploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
T0*z
_classp
nlloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::

Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1Switch§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switchloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0*z
_classp
nlloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::
č
Źloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperationDenseToDenseSetOperation loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat*
T0*<
_output_shapes*
(:˙˙˙˙˙˙˙˙˙:˙˙˙˙˙˙˙˙˙:*
set_operationa-b
˙
¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dimsSizeŽloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1*
_output_shapes
: *
T0
ç
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/xConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
dtype0*
_output_shapes
: *
value	B : 
ý
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dimsEqualloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims*
_output_shapes
: *
T0

loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1Switchloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0
*Ľ
_class
loc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
_output_shapes
: : 
â
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/MergeMergeloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims*
T0
*
N*
_output_shapes
: : 
Ł
nloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/MergeMergeloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Mergesloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1*
T0
*
N*
_output_shapes
: : 
Ç
_loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/ConstConst*8
value/B- B'weights can not be broadcast to values.*
dtype0*
_output_shapes
: 
°
aloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_1Const*
valueB Bweights.shape=*
dtype0*
_output_shapes
: 
ş
aloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_2Const*)
value B Bdense_1_sample_weights:0*
dtype0*
_output_shapes
: 
Ż
aloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_3Const*
valueB Bvalues.shape=*
dtype0*
_output_shapes
: 
Í
aloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_4Const*<
value3B1 B+loss/dense_1_loss/mean_squared_error/Mean:0*
dtype0*
_output_shapes
: 
Ź
aloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_5Const*
valueB B
is_scalar=*
dtype0*
_output_shapes
: 
ů
lloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/SwitchSwitchnloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Mergenloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge*
_output_shapes
: : *
T0


nloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_tIdentitynloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Switch:1*
T0
*
_output_shapes
: 

nloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_fIdentitylloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Switch*
T0
*
_output_shapes
: 

mloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_idIdentitynloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge*
_output_shapes
: *
T0

ă
jloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/NoOpNoOpo^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t

xloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependencyIdentitynloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_tk^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/NoOp*
T0
*
_classw
usloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t*
_output_shapes
: 
Ě
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*8
value/B- B'weights can not be broadcast to values.*
dtype0*
_output_shapes
: 
ł
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
valueB Bweights.shape=*
dtype0*
_output_shapes
: 
˝
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*)
value B Bdense_1_sample_weights:0*
dtype0*
_output_shapes
: 
˛
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
valueB Bvalues.shape=*
dtype0*
_output_shapes
: 
Đ
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
dtype0*
_output_shapes
: *<
value3B1 B+loss/dense_1_loss/mean_squared_error/Mean:0
Ż
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
dtype0*
_output_shapes
: *
valueB B
is_scalar=
˘

lloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/AssertAssertsloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switchsloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3*
T
2	


sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/SwitchSwitchnloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Mergemloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0
*
_classw
usloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge*
_output_shapes
: : 
ţ
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1Switchgloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapemloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0*z
_classp
nlloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::
ü
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2Switchfloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapemloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0*y
_classo
mkloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::
î
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3Switchcloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarmloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
_output_shapes
: : *
T0
*v
_classl
jhloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar

zloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1Identitynloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_fm^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert*
_output_shapes
: *
T0
*
_classw
usloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f

kloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/MergeMergezloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1xloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency*
N*
_output_shapes
: : *
T0


Tloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ShapeShape)loss/dense_1_loss/mean_squared_error/Meanl^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Merge*
_output_shapes
:*
T0

Tloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ConstConstl^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Merge*
valueB
 *  ?*
dtype0*
_output_shapes
: 
°
Nloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_likeFillTloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ShapeTloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/Const*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
á
Dloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weightsMuldense_1_sample_weightsNloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
Ü
6loss/dense_1_loss/mean_squared_error/weighted_loss/MulMul)loss/dense_1_loss/mean_squared_error/MeanDloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
a
loss/dense_1_loss/ConstConst*
valueB: *
dtype0*
_output_shapes
:

loss/dense_1_loss/SumSum6loss/dense_1_loss/mean_squared_error/weighted_loss/Mulloss/dense_1_loss/Const*
_output_shapes
: *
T0

loss/dense_1_loss/num_elementsSize6loss/dense_1_loss/mean_squared_error/weighted_loss/Mul*
T0*
_output_shapes
: 
{
#loss/dense_1_loss/num_elements/CastCastloss/dense_1_loss/num_elements*

SrcT0*
_output_shapes
: *

DstT0
\
loss/dense_1_loss/mul/xConst*
_output_shapes
: *
valueB
 *  ?*
dtype0
{
loss/dense_1_loss/mulMulloss/dense_1_loss/mul/x#loss/dense_1_loss/num_elements/Cast*
T0*
_output_shapes
: 
\
loss/dense_1_loss/Const_1Const*
valueB *
dtype0*
_output_shapes
: 
q
loss/dense_1_loss/Sum_1Sumloss/dense_1_loss/Sumloss/dense_1_loss/Const_1*
_output_shapes
: *
T0
t
loss/dense_1_loss/valueDivNoNanloss/dense_1_loss/Sum_1loss/dense_1_loss/mul*
_output_shapes
: *
T0
O

loss/mul/xConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 
U
loss/mulMul
loss/mul/xloss/dense_1_loss/value*
_output_shapes
: *
T0
q
iter/Initializer/zerosConst*
_output_shapes
: *
value	B	 R *
_class
	loc:@iter*
dtype0	

iterVarHandleOp"/device:CPU:0*
_class
	loc:@iter*
dtype0	*
_output_shapes
: *
shape: *
shared_nameiter
h
%iter/IsInitialized/VarIsInitializedOpVarIsInitializedOpiter"/device:CPU:0*
_output_shapes
: 
r
iter/AssignAssignVariableOpiteriter/Initializer/zeros"/device:CPU:0*
_class
	loc:@iter*
dtype0	
}
iter/Read/ReadVariableOpReadVariableOpiter"/device:CPU:0*
_class
	loc:@iter*
dtype0	*
_output_shapes
: 

'learning_rate/Initializer/initial_valueConst*
valueB
 *o:* 
_class
loc:@learning_rate*
dtype0*
_output_shapes
: 

learning_rateVarHandleOp*
dtype0*
_output_shapes
: *
shape: *
shared_namelearning_rate* 
_class
loc:@learning_rate
k
.learning_rate/IsInitialized/VarIsInitializedOpVarIsInitializedOplearning_rate*
_output_shapes
: 

learning_rate/AssignAssignVariableOplearning_rate'learning_rate/Initializer/initial_value*
dtype0* 
_class
loc:@learning_rate

!learning_rate/Read/ReadVariableOpReadVariableOplearning_rate* 
_class
loc:@learning_rate*
dtype0*
_output_shapes
: 
~
decay/Initializer/initial_valueConst*
valueB
 *    *
_class

loc:@decay*
dtype0*
_output_shapes
: 
x
decayVarHandleOp*
shared_namedecay*
_class

loc:@decay*
dtype0*
_output_shapes
: *
shape: 
[
&decay/IsInitialized/VarIsInitializedOpVarIsInitializedOpdecay*
_output_shapes
: 
o
decay/AssignAssignVariableOpdecaydecay/Initializer/initial_value*
_class

loc:@decay*
dtype0
q
decay/Read/ReadVariableOpReadVariableOpdecay*
_class

loc:@decay*
dtype0*
_output_shapes
: 

 beta_1/Initializer/initial_valueConst*
valueB
 *fff?*
_class
loc:@beta_1*
dtype0*
_output_shapes
: 
{
beta_1VarHandleOp*
_class
loc:@beta_1*
dtype0*
_output_shapes
: *
shape: *
shared_namebeta_1
]
'beta_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpbeta_1*
_output_shapes
: 
s
beta_1/AssignAssignVariableOpbeta_1 beta_1/Initializer/initial_value*
_class
loc:@beta_1*
dtype0
t
beta_1/Read/ReadVariableOpReadVariableOpbeta_1*
_class
loc:@beta_1*
dtype0*
_output_shapes
: 

 beta_2/Initializer/initial_valueConst*
valueB
 *wž?*
_class
loc:@beta_2*
dtype0*
_output_shapes
: 
{
beta_2VarHandleOp*
_class
loc:@beta_2*
dtype0*
_output_shapes
: *
shape: *
shared_namebeta_2
]
'beta_2/IsInitialized/VarIsInitializedOpVarIsInitializedOpbeta_2*
_output_shapes
: 
s
beta_2/AssignAssignVariableOpbeta_2 beta_2/Initializer/initial_value*
_class
loc:@beta_2*
dtype0
t
beta_2/Read/ReadVariableOpReadVariableOpbeta_2*
_class
loc:@beta_2*
dtype0*
_output_shapes
: 

!epsilon/Initializer/initial_valueConst*
valueB
 *żÖ3*
_class
loc:@epsilon*
dtype0*
_output_shapes
: 
~
epsilonVarHandleOp*
_output_shapes
: *
shape: *
shared_name	epsilon*
_class
loc:@epsilon*
dtype0
_
(epsilon/IsInitialized/VarIsInitializedOpVarIsInitializedOpepsilon*
_output_shapes
: 
w
epsilon/AssignAssignVariableOpepsilon!epsilon/Initializer/initial_value*
_class
loc:@epsilon*
dtype0
w
epsilon/Read/ReadVariableOpReadVariableOpepsilon*
_class
loc:@epsilon*
dtype0*
_output_shapes
: 
`
training/Adam/gradients/ShapeConst*
valueB *
dtype0*
_output_shapes
: 
f
!training/Adam/gradients/grad_ys_0Const*
valueB
 *  ?*
dtype0*
_output_shapes
: 

training/Adam/gradients/FillFilltraining/Adam/gradients/Shape!training/Adam/gradients/grad_ys_0*
T0*
_output_shapes
: 

)training/Adam/gradients/loss/mul_grad/MulMultraining/Adam/gradients/Fillloss/dense_1_loss/value*
T0*
_output_shapes
: 
}
+training/Adam/gradients/loss/mul_grad/Mul_1Multraining/Adam/gradients/Fill
loss/mul/x*
_output_shapes
: *
T0
}
:training/Adam/gradients/loss/dense_1_loss/value_grad/ShapeConst*
dtype0*
_output_shapes
: *
valueB 

<training/Adam/gradients/loss/dense_1_loss/value_grad/Shape_1Const*
valueB *
dtype0*
_output_shapes
: 

Jtraining/Adam/gradients/loss/dense_1_loss/value_grad/BroadcastGradientArgsBroadcastGradientArgs:training/Adam/gradients/loss/dense_1_loss/value_grad/Shape<training/Adam/gradients/loss/dense_1_loss/value_grad/Shape_1*2
_output_shapes 
:˙˙˙˙˙˙˙˙˙:˙˙˙˙˙˙˙˙˙
°
?training/Adam/gradients/loss/dense_1_loss/value_grad/div_no_nanDivNoNan+training/Adam/gradients/loss/mul_grad/Mul_1loss/dense_1_loss/mul*
T0*
_output_shapes
: 
í
8training/Adam/gradients/loss/dense_1_loss/value_grad/SumSum?training/Adam/gradients/loss/dense_1_loss/value_grad/div_no_nanJtraining/Adam/gradients/loss/dense_1_loss/value_grad/BroadcastGradientArgs*
_output_shapes
: *
T0
Ţ
<training/Adam/gradients/loss/dense_1_loss/value_grad/ReshapeReshape8training/Adam/gradients/loss/dense_1_loss/value_grad/Sum:training/Adam/gradients/loss/dense_1_loss/value_grad/Shape*
T0*
_output_shapes
: 
y
8training/Adam/gradients/loss/dense_1_loss/value_grad/NegNegloss/dense_1_loss/Sum_1*
_output_shapes
: *
T0
ż
Atraining/Adam/gradients/loss/dense_1_loss/value_grad/div_no_nan_1DivNoNan8training/Adam/gradients/loss/dense_1_loss/value_grad/Negloss/dense_1_loss/mul*
_output_shapes
: *
T0
Č
Atraining/Adam/gradients/loss/dense_1_loss/value_grad/div_no_nan_2DivNoNanAtraining/Adam/gradients/loss/dense_1_loss/value_grad/div_no_nan_1loss/dense_1_loss/mul*
T0*
_output_shapes
: 
Đ
8training/Adam/gradients/loss/dense_1_loss/value_grad/mulMul+training/Adam/gradients/loss/mul_grad/Mul_1Atraining/Adam/gradients/loss/dense_1_loss/value_grad/div_no_nan_2*
_output_shapes
: *
T0
ę
:training/Adam/gradients/loss/dense_1_loss/value_grad/Sum_1Sum8training/Adam/gradients/loss/dense_1_loss/value_grad/mulLtraining/Adam/gradients/loss/dense_1_loss/value_grad/BroadcastGradientArgs:1*
_output_shapes
: *
T0
ä
>training/Adam/gradients/loss/dense_1_loss/value_grad/Reshape_1Reshape:training/Adam/gradients/loss/dense_1_loss/value_grad/Sum_1<training/Adam/gradients/loss/dense_1_loss/value_grad/Shape_1*
T0*
_output_shapes
: 

Btraining/Adam/gradients/loss/dense_1_loss/Sum_1_grad/Reshape/shapeConst*
valueB *
dtype0*
_output_shapes
: 
ę
<training/Adam/gradients/loss/dense_1_loss/Sum_1_grad/ReshapeReshape<training/Adam/gradients/loss/dense_1_loss/value_grad/ReshapeBtraining/Adam/gradients/loss/dense_1_loss/Sum_1_grad/Reshape/shape*
T0*
_output_shapes
: 
}
:training/Adam/gradients/loss/dense_1_loss/Sum_1_grad/ConstConst*
valueB *
dtype0*
_output_shapes
: 
Ü
9training/Adam/gradients/loss/dense_1_loss/Sum_1_grad/TileTile<training/Adam/gradients/loss/dense_1_loss/Sum_1_grad/Reshape:training/Adam/gradients/loss/dense_1_loss/Sum_1_grad/Const*
T0*
_output_shapes
: 

@training/Adam/gradients/loss/dense_1_loss/Sum_grad/Reshape/shapeConst*
_output_shapes
:*
valueB:*
dtype0
ç
:training/Adam/gradients/loss/dense_1_loss/Sum_grad/ReshapeReshape9training/Adam/gradients/loss/dense_1_loss/Sum_1_grad/Tile@training/Adam/gradients/loss/dense_1_loss/Sum_grad/Reshape/shape*
T0*
_output_shapes
:

8training/Adam/gradients/loss/dense_1_loss/Sum_grad/ShapeShape6loss/dense_1_loss/mean_squared_error/weighted_loss/Mul*
T0*
_output_shapes
:
ă
7training/Adam/gradients/loss/dense_1_loss/Sum_grad/TileTile:training/Adam/gradients/loss/dense_1_loss/Sum_grad/Reshape8training/Adam/gradients/loss/dense_1_loss/Sum_grad/Shape*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
˛
Ytraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/ShapeShape)loss/dense_1_loss/mean_squared_error/Mean*
T0*
_output_shapes
:
Ď
[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Shape_1ShapeDloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights*
T0*
_output_shapes
:
î
itraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/BroadcastGradientArgsBroadcastGradientArgsYtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Shape[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Shape_1*2
_output_shapes 
:˙˙˙˙˙˙˙˙˙:˙˙˙˙˙˙˙˙˙

Wtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/MulMul7training/Adam/gradients/loss/dense_1_loss/Sum_grad/TileDloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
Ĺ
Wtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/SumSumWtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Mulitraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/BroadcastGradientArgs*
_output_shapes
:*
T0
Č
[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/ReshapeReshapeWtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/SumYtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Shape*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
ň
Ytraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Mul_1Mul)loss/dense_1_loss/mean_squared_error/Mean7training/Adam/gradients/loss/dense_1_loss/Sum_grad/Tile*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
Ë
Ytraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Sum_1SumYtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Mul_1ktraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/BroadcastGradientArgs:1*
T0*
_output_shapes
:
Î
]training/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Reshape_1ReshapeYtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Sum_1[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/Shape_1*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
˛
Ltraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/ShapeShape6loss/dense_1_loss/mean_squared_error/SquaredDifference*
_output_shapes
:*
T0
î
Ktraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/SizeConst*
value	B :*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
dtype0*
_output_shapes
: 
Ý
Jtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/addAdd;loss/dense_1_loss/mean_squared_error/Mean/reduction_indicesKtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Size*
T0*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
_output_shapes
: 
ń
Jtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/modFloorModJtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/addKtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Size*
T0*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
_output_shapes
: 
ň
Ntraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape_1Const*
valueB *_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
dtype0*
_output_shapes
: 
ő
Rtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/range/startConst*
dtype0*
_output_shapes
: *
value	B : *_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape
ő
Rtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/range/deltaConst*
value	B :*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
dtype0*
_output_shapes
: 
Ç
Ltraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/rangeRangeRtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/range/startKtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/SizeRtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/range/delta*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
_output_shapes
:
ô
Qtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Fill/valueConst*
value	B :*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
dtype0*
_output_shapes
: 
ř
Ktraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/FillFillNtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape_1Qtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Fill/value*
T0*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
_output_shapes
: 
Š
Ttraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/DynamicStitchDynamicStitchLtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/rangeJtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/modLtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/ShapeKtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Fill*
T0*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
N*
_output_shapes
:
ó
Ptraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Maximum/yConst*
value	B :*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
dtype0*
_output_shapes
: 

Ntraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/MaximumMaximumTtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/DynamicStitchPtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Maximum/y*
_output_shapes
:*
T0*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape
˙
Otraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/floordivFloorDivLtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/ShapeNtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Maximum*
T0*_
_classU
SQloc:@training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape*
_output_shapes
:
Ç
Ntraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/ReshapeReshape[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/weighted_loss/Mul_grad/ReshapeTtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/DynamicStitch*
T0*0
_output_shapes
:˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙
Ż
Ktraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/TileTileNtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/ReshapeOtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/floordiv*
T0*0
_output_shapes
:˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙
´
Ntraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape_2Shape6loss/dense_1_loss/mean_squared_error/SquaredDifference*
T0*
_output_shapes
:
§
Ntraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape_3Shape)loss/dense_1_loss/mean_squared_error/Mean*
T0*
_output_shapes
:

Ltraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/ConstConst*
valueB: *
dtype0*
_output_shapes
:

Ktraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/ProdProdNtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape_2Ltraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Const*
T0*
_output_shapes
: 

Ntraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Const_1Const*
valueB: *
dtype0*
_output_shapes
:

Mtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Prod_1ProdNtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Shape_3Ntraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Const_1*
T0*
_output_shapes
: 

Rtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Maximum_1/yConst*
_output_shapes
: *
value	B :*
dtype0

Ptraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Maximum_1MaximumMtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Prod_1Rtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Maximum_1/y*
_output_shapes
: *
T0

Qtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/floordiv_1FloorDivKtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/ProdPtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Maximum_1*
T0*
_output_shapes
: 
Ö
Ktraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/CastCastQtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/floordiv_1*

SrcT0*
_output_shapes
: *

DstT0
Ľ
Ntraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/truedivRealDivKtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/TileKtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/Cast*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0

Ytraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/ShapeShapedense_1/BiasAdd*
T0*
_output_shapes
:

[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Shape_1Shapedense_1_target*
T0*
_output_shapes
:
î
itraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/BroadcastGradientArgsBroadcastGradientArgsYtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Shape[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Shape_1*2
_output_shapes 
:˙˙˙˙˙˙˙˙˙:˙˙˙˙˙˙˙˙˙
đ
Ztraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/scalarConstO^training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/truediv*
valueB
 *   @*
dtype0*
_output_shapes
: 
ź
Wtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/MulMulZtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/scalarNtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/truediv*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙

Wtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/subSubdense_1/BiasAdddense_1_targetO^training/Adam/gradients/loss/dense_1_loss/mean_squared_error/Mean_grad/truediv*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
Ä
Ytraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/mul_1MulWtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/MulWtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/sub*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
Ç
Wtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/SumSumYtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/mul_1itraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/BroadcastGradientArgs*
_output_shapes
:*
T0
Ě
[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/ReshapeReshapeWtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/SumYtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Shape*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
Ë
Ytraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Sum_1SumYtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/mul_1ktraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/BroadcastGradientArgs:1*
T0*
_output_shapes
:
Ű
]training/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Reshape_1ReshapeYtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Sum_1[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Shape_1*
T0*0
_output_shapes
:˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙
ř
Wtraining/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/NegNeg]training/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Reshape_1*0
_output_shapes
:˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙*
T0
É
8training/Adam/gradients/dense_1/BiasAdd_grad/BiasAddGradBiasAddGrad[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Reshape*
_output_shapes
:*
T0
ý
2training/Adam/gradients/dense_1/MatMul_grad/MatMulMatMul[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Reshapedense_1/MatMul/ReadVariableOp*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
transpose_b(
ă
4training/Adam/gradients/dense_1/MatMul_grad/MatMul_1MatMul
dense/Relu[training/Adam/gradients/loss/dense_1_loss/mean_squared_error/SquaredDifference_grad/Reshape*
T0*
_output_shapes

:*
transpose_a(
Ž
0training/Adam/gradients/dense/Relu_grad/ReluGradReluGrad2training/Adam/gradients/dense_1/MatMul_grad/MatMul
dense/Relu*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙

6training/Adam/gradients/dense/BiasAdd_grad/BiasAddGradBiasAddGrad0training/Adam/gradients/dense/Relu_grad/ReluGrad*
_output_shapes
:*
T0
Î
0training/Adam/gradients/dense/MatMul_grad/MatMulMatMul0training/Adam/gradients/dense/Relu_grad/ReluGraddense/MatMul/ReadVariableOp*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
*
transpose_b(*
T0
ˇ
2training/Adam/gradients/dense/MatMul_grad/MatMul_1MatMuldense_input0training/Adam/gradients/dense/Relu_grad/ReluGrad*
T0*
_output_shapes

:
*
transpose_a(
´
.training/Adam/dense/kernel/m/Initializer/zerosConst*
dtype0*
_output_shapes

:
*
valueB
*    */
_class%
#!loc:@training/Adam/dense/kernel/m
Ĺ
training/Adam/dense/kernel/mVarHandleOp*
shape
:
*-
shared_nametraining/Adam/dense/kernel/m*/
_class%
#!loc:@training/Adam/dense/kernel/m*
dtype0*
_output_shapes
: 

=training/Adam/dense/kernel/m/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense/kernel/m*
_output_shapes
: 
Ă
#training/Adam/dense/kernel/m/AssignAssignVariableOptraining/Adam/dense/kernel/m.training/Adam/dense/kernel/m/Initializer/zeros*/
_class%
#!loc:@training/Adam/dense/kernel/m*
dtype0
ž
0training/Adam/dense/kernel/m/Read/ReadVariableOpReadVariableOptraining/Adam/dense/kernel/m*/
_class%
#!loc:@training/Adam/dense/kernel/m*
dtype0*
_output_shapes

:

¨
,training/Adam/dense/bias/m/Initializer/zerosConst*
dtype0*
_output_shapes
:*
valueB*    *-
_class#
!loc:@training/Adam/dense/bias/m
ť
training/Adam/dense/bias/mVarHandleOp*-
_class#
!loc:@training/Adam/dense/bias/m*
dtype0*
_output_shapes
: *
shape:*+
shared_nametraining/Adam/dense/bias/m

;training/Adam/dense/bias/m/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense/bias/m*
_output_shapes
: 
ť
!training/Adam/dense/bias/m/AssignAssignVariableOptraining/Adam/dense/bias/m,training/Adam/dense/bias/m/Initializer/zeros*-
_class#
!loc:@training/Adam/dense/bias/m*
dtype0
´
.training/Adam/dense/bias/m/Read/ReadVariableOpReadVariableOptraining/Adam/dense/bias/m*-
_class#
!loc:@training/Adam/dense/bias/m*
dtype0*
_output_shapes
:
¸
0training/Adam/dense_1/kernel/m/Initializer/zerosConst*
_output_shapes

:*
valueB*    *1
_class'
%#loc:@training/Adam/dense_1/kernel/m*
dtype0
Ë
training/Adam/dense_1/kernel/mVarHandleOp*
shape
:*/
shared_name training/Adam/dense_1/kernel/m*1
_class'
%#loc:@training/Adam/dense_1/kernel/m*
dtype0*
_output_shapes
: 

?training/Adam/dense_1/kernel/m/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense_1/kernel/m*
_output_shapes
: 
Ë
%training/Adam/dense_1/kernel/m/AssignAssignVariableOptraining/Adam/dense_1/kernel/m0training/Adam/dense_1/kernel/m/Initializer/zeros*1
_class'
%#loc:@training/Adam/dense_1/kernel/m*
dtype0
Ä
2training/Adam/dense_1/kernel/m/Read/ReadVariableOpReadVariableOptraining/Adam/dense_1/kernel/m*1
_class'
%#loc:@training/Adam/dense_1/kernel/m*
dtype0*
_output_shapes

:
Ź
.training/Adam/dense_1/bias/m/Initializer/zerosConst*
valueB*    */
_class%
#!loc:@training/Adam/dense_1/bias/m*
dtype0*
_output_shapes
:
Á
training/Adam/dense_1/bias/mVarHandleOp*
shape:*-
shared_nametraining/Adam/dense_1/bias/m*/
_class%
#!loc:@training/Adam/dense_1/bias/m*
dtype0*
_output_shapes
: 

=training/Adam/dense_1/bias/m/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense_1/bias/m*
_output_shapes
: 
Ă
#training/Adam/dense_1/bias/m/AssignAssignVariableOptraining/Adam/dense_1/bias/m.training/Adam/dense_1/bias/m/Initializer/zeros*/
_class%
#!loc:@training/Adam/dense_1/bias/m*
dtype0
ş
0training/Adam/dense_1/bias/m/Read/ReadVariableOpReadVariableOptraining/Adam/dense_1/bias/m*/
_class%
#!loc:@training/Adam/dense_1/bias/m*
dtype0*
_output_shapes
:
´
.training/Adam/dense/kernel/v/Initializer/zerosConst*
valueB
*    */
_class%
#!loc:@training/Adam/dense/kernel/v*
dtype0*
_output_shapes

:

Ĺ
training/Adam/dense/kernel/vVarHandleOp*/
_class%
#!loc:@training/Adam/dense/kernel/v*
dtype0*
_output_shapes
: *
shape
:
*-
shared_nametraining/Adam/dense/kernel/v

=training/Adam/dense/kernel/v/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense/kernel/v*
_output_shapes
: 
Ă
#training/Adam/dense/kernel/v/AssignAssignVariableOptraining/Adam/dense/kernel/v.training/Adam/dense/kernel/v/Initializer/zeros*/
_class%
#!loc:@training/Adam/dense/kernel/v*
dtype0
ž
0training/Adam/dense/kernel/v/Read/ReadVariableOpReadVariableOptraining/Adam/dense/kernel/v*/
_class%
#!loc:@training/Adam/dense/kernel/v*
dtype0*
_output_shapes

:

¨
,training/Adam/dense/bias/v/Initializer/zerosConst*
valueB*    *-
_class#
!loc:@training/Adam/dense/bias/v*
dtype0*
_output_shapes
:
ť
training/Adam/dense/bias/vVarHandleOp*
shape:*+
shared_nametraining/Adam/dense/bias/v*-
_class#
!loc:@training/Adam/dense/bias/v*
dtype0*
_output_shapes
: 

;training/Adam/dense/bias/v/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense/bias/v*
_output_shapes
: 
ť
!training/Adam/dense/bias/v/AssignAssignVariableOptraining/Adam/dense/bias/v,training/Adam/dense/bias/v/Initializer/zeros*
dtype0*-
_class#
!loc:@training/Adam/dense/bias/v
´
.training/Adam/dense/bias/v/Read/ReadVariableOpReadVariableOptraining/Adam/dense/bias/v*
_output_shapes
:*-
_class#
!loc:@training/Adam/dense/bias/v*
dtype0
¸
0training/Adam/dense_1/kernel/v/Initializer/zerosConst*
valueB*    *1
_class'
%#loc:@training/Adam/dense_1/kernel/v*
dtype0*
_output_shapes

:
Ë
training/Adam/dense_1/kernel/vVarHandleOp*1
_class'
%#loc:@training/Adam/dense_1/kernel/v*
dtype0*
_output_shapes
: *
shape
:*/
shared_name training/Adam/dense_1/kernel/v

?training/Adam/dense_1/kernel/v/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense_1/kernel/v*
_output_shapes
: 
Ë
%training/Adam/dense_1/kernel/v/AssignAssignVariableOptraining/Adam/dense_1/kernel/v0training/Adam/dense_1/kernel/v/Initializer/zeros*1
_class'
%#loc:@training/Adam/dense_1/kernel/v*
dtype0
Ä
2training/Adam/dense_1/kernel/v/Read/ReadVariableOpReadVariableOptraining/Adam/dense_1/kernel/v*1
_class'
%#loc:@training/Adam/dense_1/kernel/v*
dtype0*
_output_shapes

:
Ź
.training/Adam/dense_1/bias/v/Initializer/zerosConst*
valueB*    */
_class%
#!loc:@training/Adam/dense_1/bias/v*
dtype0*
_output_shapes
:
Á
training/Adam/dense_1/bias/vVarHandleOp*-
shared_nametraining/Adam/dense_1/bias/v*/
_class%
#!loc:@training/Adam/dense_1/bias/v*
dtype0*
_output_shapes
: *
shape:

=training/Adam/dense_1/bias/v/IsInitialized/VarIsInitializedOpVarIsInitializedOptraining/Adam/dense_1/bias/v*
_output_shapes
: 
Ă
#training/Adam/dense_1/bias/v/AssignAssignVariableOptraining/Adam/dense_1/bias/v.training/Adam/dense_1/bias/v/Initializer/zeros*/
_class%
#!loc:@training/Adam/dense_1/bias/v*
dtype0
ş
0training/Adam/dense_1/bias/v/Read/ReadVariableOpReadVariableOptraining/Adam/dense_1/bias/v*/
_class%
#!loc:@training/Adam/dense_1/bias/v*
dtype0*
_output_shapes
:

5training/Adam/Adam/update_dense/kernel/ReadVariableOpReadVariableOpiter"/device:CPU:0*
dtype0	*
_output_shapes
: 

,training/Adam/Adam/update_dense/kernel/add/yConst*
dtype0	*
_output_shapes
: *
value	B	 R*
_class
loc:@dense/kernel
Ř
*training/Adam/Adam/update_dense/kernel/addAdd5training/Adam/Adam/update_dense/kernel/ReadVariableOp,training/Adam/Adam/update_dense/kernel/add/y*
T0	*
_class
loc:@dense/kernel*
_output_shapes
: 
°
+training/Adam/Adam/update_dense/kernel/CastCast*training/Adam/Adam/update_dense/kernel/add*

SrcT0	*
_class
loc:@dense/kernel*
_output_shapes
: *

DstT0
x
9training/Adam/Adam/update_dense/kernel/Pow/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: 
Ű
*training/Adam/Adam/update_dense/kernel/PowPow9training/Adam/Adam/update_dense/kernel/Pow/ReadVariableOp+training/Adam/Adam/update_dense/kernel/Cast*
T0*
_class
loc:@dense/kernel*
_output_shapes
: 
z
;training/Adam/Adam/update_dense/kernel/Pow_1/ReadVariableOpReadVariableOpbeta_2*
_output_shapes
: *
dtype0
ß
,training/Adam/Adam/update_dense/kernel/Pow_1Pow;training/Adam/Adam/update_dense/kernel/Pow_1/ReadVariableOp+training/Adam/Adam/update_dense/kernel/Cast*
T0*
_class
loc:@dense/kernel*
_output_shapes
: 

Gtraining/Adam/Adam/update_dense/kernel/ResourceApplyAdam/ReadVariableOpReadVariableOplearning_rate*
dtype0*
_output_shapes
: 

Itraining/Adam/Adam/update_dense/kernel/ResourceApplyAdam/ReadVariableOp_1ReadVariableOpbeta_1*
dtype0*
_output_shapes
: 

Itraining/Adam/Adam/update_dense/kernel/ResourceApplyAdam/ReadVariableOp_2ReadVariableOpbeta_2*
_output_shapes
: *
dtype0

Itraining/Adam/Adam/update_dense/kernel/ResourceApplyAdam/ReadVariableOp_3ReadVariableOpepsilon*
dtype0*
_output_shapes
: 

8training/Adam/Adam/update_dense/kernel/ResourceApplyAdamResourceApplyAdamdense/kerneltraining/Adam/dense/kernel/mtraining/Adam/dense/kernel/v*training/Adam/Adam/update_dense/kernel/Pow,training/Adam/Adam/update_dense/kernel/Pow_1Gtraining/Adam/Adam/update_dense/kernel/ResourceApplyAdam/ReadVariableOpItraining/Adam/Adam/update_dense/kernel/ResourceApplyAdam/ReadVariableOp_1Itraining/Adam/Adam/update_dense/kernel/ResourceApplyAdam/ReadVariableOp_2Itraining/Adam/Adam/update_dense/kernel/ResourceApplyAdam/ReadVariableOp_32training/Adam/gradients/dense/MatMul_grad/MatMul_1*
use_locking(*
T0*
_class
loc:@dense/kernel

3training/Adam/Adam/update_dense/bias/ReadVariableOpReadVariableOpiter"/device:CPU:0*
dtype0	*
_output_shapes
: 

*training/Adam/Adam/update_dense/bias/add/yConst*
value	B	 R*
_class
loc:@dense/bias*
dtype0	*
_output_shapes
: 
Đ
(training/Adam/Adam/update_dense/bias/addAdd3training/Adam/Adam/update_dense/bias/ReadVariableOp*training/Adam/Adam/update_dense/bias/add/y*
_class
loc:@dense/bias*
_output_shapes
: *
T0	
Ş
)training/Adam/Adam/update_dense/bias/CastCast(training/Adam/Adam/update_dense/bias/add*

SrcT0	*
_class
loc:@dense/bias*
_output_shapes
: *

DstT0
v
7training/Adam/Adam/update_dense/bias/Pow/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: 
Ó
(training/Adam/Adam/update_dense/bias/PowPow7training/Adam/Adam/update_dense/bias/Pow/ReadVariableOp)training/Adam/Adam/update_dense/bias/Cast*
T0*
_class
loc:@dense/bias*
_output_shapes
: 
x
9training/Adam/Adam/update_dense/bias/Pow_1/ReadVariableOpReadVariableOpbeta_2*
dtype0*
_output_shapes
: 
×
*training/Adam/Adam/update_dense/bias/Pow_1Pow9training/Adam/Adam/update_dense/bias/Pow_1/ReadVariableOp)training/Adam/Adam/update_dense/bias/Cast*
T0*
_class
loc:@dense/bias*
_output_shapes
: 

Etraining/Adam/Adam/update_dense/bias/ResourceApplyAdam/ReadVariableOpReadVariableOplearning_rate*
dtype0*
_output_shapes
: 

Gtraining/Adam/Adam/update_dense/bias/ResourceApplyAdam/ReadVariableOp_1ReadVariableOpbeta_1*
dtype0*
_output_shapes
: 

Gtraining/Adam/Adam/update_dense/bias/ResourceApplyAdam/ReadVariableOp_2ReadVariableOpbeta_2*
_output_shapes
: *
dtype0

Gtraining/Adam/Adam/update_dense/bias/ResourceApplyAdam/ReadVariableOp_3ReadVariableOpepsilon*
dtype0*
_output_shapes
: 
ú
6training/Adam/Adam/update_dense/bias/ResourceApplyAdamResourceApplyAdam
dense/biastraining/Adam/dense/bias/mtraining/Adam/dense/bias/v(training/Adam/Adam/update_dense/bias/Pow*training/Adam/Adam/update_dense/bias/Pow_1Etraining/Adam/Adam/update_dense/bias/ResourceApplyAdam/ReadVariableOpGtraining/Adam/Adam/update_dense/bias/ResourceApplyAdam/ReadVariableOp_1Gtraining/Adam/Adam/update_dense/bias/ResourceApplyAdam/ReadVariableOp_2Gtraining/Adam/Adam/update_dense/bias/ResourceApplyAdam/ReadVariableOp_36training/Adam/gradients/dense/BiasAdd_grad/BiasAddGrad*
use_locking(*
T0*
_class
loc:@dense/bias

7training/Adam/Adam/update_dense_1/kernel/ReadVariableOpReadVariableOpiter"/device:CPU:0*
dtype0	*
_output_shapes
: 

.training/Adam/Adam/update_dense_1/kernel/add/yConst*
dtype0	*
_output_shapes
: *
value	B	 R*!
_class
loc:@dense_1/kernel
ŕ
,training/Adam/Adam/update_dense_1/kernel/addAdd7training/Adam/Adam/update_dense_1/kernel/ReadVariableOp.training/Adam/Adam/update_dense_1/kernel/add/y*
T0	*!
_class
loc:@dense_1/kernel*
_output_shapes
: 
ś
-training/Adam/Adam/update_dense_1/kernel/CastCast,training/Adam/Adam/update_dense_1/kernel/add*

SrcT0	*!
_class
loc:@dense_1/kernel*
_output_shapes
: *

DstT0
z
;training/Adam/Adam/update_dense_1/kernel/Pow/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: 
ă
,training/Adam/Adam/update_dense_1/kernel/PowPow;training/Adam/Adam/update_dense_1/kernel/Pow/ReadVariableOp-training/Adam/Adam/update_dense_1/kernel/Cast*
_output_shapes
: *
T0*!
_class
loc:@dense_1/kernel
|
=training/Adam/Adam/update_dense_1/kernel/Pow_1/ReadVariableOpReadVariableOpbeta_2*
dtype0*
_output_shapes
: 
ç
.training/Adam/Adam/update_dense_1/kernel/Pow_1Pow=training/Adam/Adam/update_dense_1/kernel/Pow_1/ReadVariableOp-training/Adam/Adam/update_dense_1/kernel/Cast*
T0*!
_class
loc:@dense_1/kernel*
_output_shapes
: 

Itraining/Adam/Adam/update_dense_1/kernel/ResourceApplyAdam/ReadVariableOpReadVariableOplearning_rate*
dtype0*
_output_shapes
: 

Ktraining/Adam/Adam/update_dense_1/kernel/ResourceApplyAdam/ReadVariableOp_1ReadVariableOpbeta_1*
dtype0*
_output_shapes
: 

Ktraining/Adam/Adam/update_dense_1/kernel/ResourceApplyAdam/ReadVariableOp_2ReadVariableOpbeta_2*
dtype0*
_output_shapes
: 

Ktraining/Adam/Adam/update_dense_1/kernel/ResourceApplyAdam/ReadVariableOp_3ReadVariableOpepsilon*
dtype0*
_output_shapes
: 
¤
:training/Adam/Adam/update_dense_1/kernel/ResourceApplyAdamResourceApplyAdamdense_1/kerneltraining/Adam/dense_1/kernel/mtraining/Adam/dense_1/kernel/v,training/Adam/Adam/update_dense_1/kernel/Pow.training/Adam/Adam/update_dense_1/kernel/Pow_1Itraining/Adam/Adam/update_dense_1/kernel/ResourceApplyAdam/ReadVariableOpKtraining/Adam/Adam/update_dense_1/kernel/ResourceApplyAdam/ReadVariableOp_1Ktraining/Adam/Adam/update_dense_1/kernel/ResourceApplyAdam/ReadVariableOp_2Ktraining/Adam/Adam/update_dense_1/kernel/ResourceApplyAdam/ReadVariableOp_34training/Adam/gradients/dense_1/MatMul_grad/MatMul_1*
use_locking(*
T0*!
_class
loc:@dense_1/kernel

5training/Adam/Adam/update_dense_1/bias/ReadVariableOpReadVariableOpiter"/device:CPU:0*
dtype0	*
_output_shapes
: 

,training/Adam/Adam/update_dense_1/bias/add/yConst*
_output_shapes
: *
value	B	 R*
_class
loc:@dense_1/bias*
dtype0	
Ř
*training/Adam/Adam/update_dense_1/bias/addAdd5training/Adam/Adam/update_dense_1/bias/ReadVariableOp,training/Adam/Adam/update_dense_1/bias/add/y*
T0	*
_class
loc:@dense_1/bias*
_output_shapes
: 
°
+training/Adam/Adam/update_dense_1/bias/CastCast*training/Adam/Adam/update_dense_1/bias/add*

SrcT0	*
_class
loc:@dense_1/bias*
_output_shapes
: *

DstT0
x
9training/Adam/Adam/update_dense_1/bias/Pow/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: 
Ű
*training/Adam/Adam/update_dense_1/bias/PowPow9training/Adam/Adam/update_dense_1/bias/Pow/ReadVariableOp+training/Adam/Adam/update_dense_1/bias/Cast*
_output_shapes
: *
T0*
_class
loc:@dense_1/bias
z
;training/Adam/Adam/update_dense_1/bias/Pow_1/ReadVariableOpReadVariableOpbeta_2*
dtype0*
_output_shapes
: 
ß
,training/Adam/Adam/update_dense_1/bias/Pow_1Pow;training/Adam/Adam/update_dense_1/bias/Pow_1/ReadVariableOp+training/Adam/Adam/update_dense_1/bias/Cast*
T0*
_class
loc:@dense_1/bias*
_output_shapes
: 

Gtraining/Adam/Adam/update_dense_1/bias/ResourceApplyAdam/ReadVariableOpReadVariableOplearning_rate*
dtype0*
_output_shapes
: 

Itraining/Adam/Adam/update_dense_1/bias/ResourceApplyAdam/ReadVariableOp_1ReadVariableOpbeta_1*
dtype0*
_output_shapes
: 

Itraining/Adam/Adam/update_dense_1/bias/ResourceApplyAdam/ReadVariableOp_2ReadVariableOpbeta_2*
dtype0*
_output_shapes
: 

Itraining/Adam/Adam/update_dense_1/bias/ResourceApplyAdam/ReadVariableOp_3ReadVariableOpepsilon*
dtype0*
_output_shapes
: 

8training/Adam/Adam/update_dense_1/bias/ResourceApplyAdamResourceApplyAdamdense_1/biastraining/Adam/dense_1/bias/mtraining/Adam/dense_1/bias/v*training/Adam/Adam/update_dense_1/bias/Pow,training/Adam/Adam/update_dense_1/bias/Pow_1Gtraining/Adam/Adam/update_dense_1/bias/ResourceApplyAdam/ReadVariableOpItraining/Adam/Adam/update_dense_1/bias/ResourceApplyAdam/ReadVariableOp_1Itraining/Adam/Adam/update_dense_1/bias/ResourceApplyAdam/ReadVariableOp_2Itraining/Adam/Adam/update_dense_1/bias/ResourceApplyAdam/ReadVariableOp_38training/Adam/gradients/dense_1/BiasAdd_grad/BiasAddGrad*
use_locking(*
T0*
_class
loc:@dense_1/bias
Ć
training/Adam/Adam/ConstConst7^training/Adam/Adam/update_dense/bias/ResourceApplyAdam9^training/Adam/Adam/update_dense/kernel/ResourceApplyAdam9^training/Adam/Adam/update_dense_1/bias/ResourceApplyAdam;^training/Adam/Adam/update_dense_1/kernel/ResourceApplyAdam*
dtype0	*
_output_shapes
: *
value	B	 R
j
&training/Adam/Adam/AssignAddVariableOpAssignAddVariableOpitertraining/Adam/Adam/Const*
dtype0	
ó
!training/Adam/Adam/ReadVariableOpReadVariableOpiter'^training/Adam/Adam/AssignAddVariableOp7^training/Adam/Adam/update_dense/bias/ResourceApplyAdam9^training/Adam/Adam/update_dense/kernel/ResourceApplyAdam9^training/Adam/Adam/update_dense_1/bias/ResourceApplyAdam;^training/Adam/Adam/update_dense_1/kernel/ResourceApplyAdam*
dtype0	*
_output_shapes
: 
i
training_1/group_depsNoOp	^loss/mul^metrics/accuracy/Mean'^training/Adam/Adam/AssignAddVariableOp
W
Const_1Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_2Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_3Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_4Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
F
VarIsInitializedOpVarIsInitializedOpiter*
_output_shapes
: 
N
VarIsInitializedOp_1VarIsInitializedOp
dense/bias*
_output_shapes
: 
Q
VarIsInitializedOp_2VarIsInitializedOplearning_rate*
_output_shapes
: 
K
VarIsInitializedOp_3VarIsInitializedOpepsilon*
_output_shapes
: 
I
VarIsInitializedOp_4VarIsInitializedOpdecay*
_output_shapes
: 
R
VarIsInitializedOp_5VarIsInitializedOpdense_1/kernel*
_output_shapes
: 
`
VarIsInitializedOp_6VarIsInitializedOptraining/Adam/dense_1/bias/m*
_output_shapes
: 
I
VarIsInitializedOp_7VarIsInitializedOptotal*
_output_shapes
: 
b
VarIsInitializedOp_8VarIsInitializedOptraining/Adam/dense_1/kernel/v*
_output_shapes
: 
`
VarIsInitializedOp_9VarIsInitializedOptraining/Adam/dense_1/bias/v*
_output_shapes
: 
a
VarIsInitializedOp_10VarIsInitializedOptraining/Adam/dense/kernel/m*
_output_shapes
: 
J
VarIsInitializedOp_11VarIsInitializedOpcount*
_output_shapes
: 
a
VarIsInitializedOp_12VarIsInitializedOptraining/Adam/dense/kernel/v*
_output_shapes
: 
_
VarIsInitializedOp_13VarIsInitializedOptraining/Adam/dense/bias/v*
_output_shapes
: 
K
VarIsInitializedOp_14VarIsInitializedOpbeta_1*
_output_shapes
: 
_
VarIsInitializedOp_15VarIsInitializedOptraining/Adam/dense/bias/m*
_output_shapes
: 
Q
VarIsInitializedOp_16VarIsInitializedOpdense_1/bias*
_output_shapes
: 
c
VarIsInitializedOp_17VarIsInitializedOptraining/Adam/dense_1/kernel/m*
_output_shapes
: 
Q
VarIsInitializedOp_18VarIsInitializedOpdense/kernel*
_output_shapes
: 
K
VarIsInitializedOp_19VarIsInitializedOpbeta_2*
_output_shapes
: 

	init/NoOpNoOp^beta_1/Assign^beta_2/Assign^count/Assign^decay/Assign^dense/bias/Assign^dense/kernel/Assign^dense_1/bias/Assign^dense_1/kernel/Assign^epsilon/Assign^learning_rate/Assign^total/Assign"^training/Adam/dense/bias/m/Assign"^training/Adam/dense/bias/v/Assign$^training/Adam/dense/kernel/m/Assign$^training/Adam/dense/kernel/v/Assign$^training/Adam/dense_1/bias/m/Assign$^training/Adam/dense_1/bias/v/Assign&^training/Adam/dense_1/kernel/m/Assign&^training/Adam/dense_1/kernel/v/Assign
0
init/NoOp_1NoOp^iter/Assign"/device:CPU:0
&
initNoOp
^init/NoOp^init/NoOp_1
W
Const_5Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
\
Const_6Const"/device:CPU:0*
valueB Bmodel*
dtype0*
_output_shapes
: 
W
Const_7Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_8Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_9Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_10Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_11Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
¤
RestoreV2/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
r
RestoreV2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

	RestoreV2	RestoreV2Const_6RestoreV2/tensor_namesRestoreV2/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
B
IdentityIdentity	RestoreV2*
_output_shapes
:*
T0
I
AssignVariableOpAssignVariableOpdense/kernelIdentity*
dtype0
¤
RestoreV2_1/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_1	RestoreV2Const_6RestoreV2_1/tensor_namesRestoreV2_1/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_1IdentityRestoreV2_1*
_output_shapes
:*
T0
K
AssignVariableOp_1AssignVariableOp
dense/bias
Identity_1*
dtype0
Ś
RestoreV2_2/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_2/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_2	RestoreV2Const_6RestoreV2_2/tensor_namesRestoreV2_2/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_2IdentityRestoreV2_2*
T0*
_output_shapes
:
O
AssignVariableOp_2AssignVariableOpdense_1/kernel
Identity_2*
dtype0
¤
RestoreV2_3/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*I
value@B>B4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE
t
RestoreV2_3/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_3	RestoreV2Const_6RestoreV2_3/tensor_namesRestoreV2_3/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_3IdentityRestoreV2_3*
T0*
_output_shapes
:
M
AssignVariableOp_3AssignVariableOpdense_1/bias
Identity_3*
dtype0
X
Const_12Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_13Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
î
SaveV2/tensor_namesConst"/device:CPU:0*
valueBB/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB(optimizer/.ATTRIBUTES/OBJECT_CONFIG_JSONB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH*
dtype0*
_output_shapes
:

SaveV2/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*C
value:B8B B B B B B B B B B B B B B B B B B B B B B B B 
ý
SaveV2SaveV2Const_13SaveV2/tensor_namesSaveV2/shape_and_slicesConst_7Const_8Const_9Const_10Const_11 dense/kernel/Read/ReadVariableOpdense/bias/Read/ReadVariableOp"dense_1/kernel/Read/ReadVariableOp dense_1/bias/Read/ReadVariableOpiter/Read/ReadVariableOp!learning_rate/Read/ReadVariableOpdecay/Read/ReadVariableOpbeta_1/Read/ReadVariableOpbeta_2/Read/ReadVariableOpepsilon/Read/ReadVariableOp0training/Adam/dense/kernel/m/Read/ReadVariableOp.training/Adam/dense/bias/m/Read/ReadVariableOp2training/Adam/dense_1/kernel/m/Read/ReadVariableOp0training/Adam/dense_1/bias/m/Read/ReadVariableOp0training/Adam/dense/kernel/v/Read/ReadVariableOp.training/Adam/dense/bias/v/Read/ReadVariableOp2training/Adam/dense_1/kernel/v/Read/ReadVariableOp0training/Adam/dense_1/bias/v/Read/ReadVariableOpConst_12"/device:CPU:0*&
dtypes
2	
Y

Identity_4IdentityConst_13^SaveV2"/device:CPU:0*
T0*
_output_shapes
: 
z
total_1/Initializer/zerosConst*
_output_shapes
: *
valueB
 *    *
_class
loc:@total_1*
dtype0
~
total_1VarHandleOp*
_output_shapes
: *
shape: *
shared_name	total_1*
_class
loc:@total_1*
dtype0
_
(total_1/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal_1*
_output_shapes
: 
o
total_1/AssignAssignVariableOptotal_1total_1/Initializer/zeros*
_class
loc:@total_1*
dtype0
w
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_class
loc:@total_1*
dtype0*
_output_shapes
: 
z
count_1/Initializer/zerosConst*
valueB
 *    *
_class
loc:@count_1*
dtype0*
_output_shapes
: 
~
count_1VarHandleOp*
_output_shapes
: *
shape: *
shared_name	count_1*
_class
loc:@count_1*
dtype0
_
(count_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount_1*
_output_shapes
: 
o
count_1/AssignAssignVariableOpcount_1count_1/Initializer/zeros*
_class
loc:@count_1*
dtype0
w
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_class
loc:@count_1*
dtype0*
_output_shapes
: 
K
Const_14Const*
dtype0*
_output_shapes
: *
valueB 
L
SumSummetrics/accuracy/MeanConst_14*
T0*
_output_shapes
: 
E
AssignAddVariableOpAssignAddVariableOptotal_1Sum*
dtype0
j
ReadVariableOpReadVariableOptotal_1^AssignAddVariableOp^Sum*
dtype0*
_output_shapes
: 
F
SizeConst*
value	B :*
dtype0*
_output_shapes
: 
B
CastCastSize*
_output_shapes
: *

DstT0*

SrcT0
^
AssignAddVariableOp_1AssignAddVariableOpcount_1Cast^AssignAddVariableOp*
dtype0
~
ReadVariableOp_1ReadVariableOpcount_1^AssignAddVariableOp^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
q
div_no_nan/ReadVariableOpReadVariableOptotal_1^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
s
div_no_nan/ReadVariableOp_1ReadVariableOpcount_1^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
o

div_no_nanDivNoNandiv_no_nan/ReadVariableOpdiv_no_nan/ReadVariableOp_1*
_output_shapes
: *
T0
C

Identity_5Identity
div_no_nan*
T0*
_output_shapes
: 
[
div_no_nan_1/ReadVariableOpReadVariableOptotal_1*
dtype0*
_output_shapes
: 
]
div_no_nan_1/ReadVariableOp_1ReadVariableOpcount_1*
dtype0*
_output_shapes
: 
u
div_no_nan_1DivNoNandiv_no_nan_1/ReadVariableOpdiv_no_nan_1/ReadVariableOp_1*
T0*
_output_shapes
: 
E

Identity_6Identitydiv_no_nan_1*
T0*
_output_shapes
: 
l
metric_op_wrapperConst^AssignAddVariableOp_1*
valueB *
dtype0*
_output_shapes
: 
Y
save/filename/inputConst*
valueB Bmodel*
dtype0*
_output_shapes
: 
n
save/filenamePlaceholderWithDefaultsave/filename/input*
_output_shapes
: *
shape: *
dtype0
e

save/ConstPlaceholderWithDefaultsave/filename*
_output_shapes
: *
shape: *
dtype0

save/Const_1Const*Í
valueĂBŔ Bš{"class_name": "Sequential", "config": {"layers": [{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "batch_input_shape": [null, 10], "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense", "trainable": true, "units": 5, "use_bias": true}}, {"class_name": "Dense", "config": {"activation": "linear", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense_1", "trainable": true, "units": 2, "use_bias": true}}], "name": "sequential"}}*
dtype0*
_output_shapes
: 
Ö
save/Const_2Const*
_output_shapes
: *
valueB B{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 10], "dtype": "float32", "name": "dense_input", "sparse": false}}*
dtype0

save/Const_3Const*
_output_shapes
: *Ř
valueÎBË BÄ{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "batch_input_shape": [null, 10], "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense", "trainable": true, "units": 5, "use_bias": true}}*
dtype0
ř
save/Const_4Const*ť
valueąBŽ B§{"class_name": "Dense", "config": {"activation": "linear", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense_1", "trainable": true, "units": 2, "use_bias": true}}*
dtype0*
_output_shapes
: 
N
save/VarIsInitializedOpVarIsInitializedOptotal_1*
_output_shapes
: 
P
save/VarIsInitializedOp_1VarIsInitializedOpcount_1*
_output_shapes
: 
3
	save/initNoOp^count_1/Assign^total_1/Assign
Ş
save/Const_5Const*
_output_shapes
: *í
valueăBŕ BŮ{"class_name": "Adam", "config": {"amsgrad": false, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "decay": 0.0, "epsilon": 1.0000000116860974e-07, "learning_rate": 0.0010000000474974513, "name": "Adam"}}*
dtype0
Ć
save/SaveV2/tensor_namesConst*ů

valueď
Bě
B/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB(optimizer/.ATTRIBUTES/OBJECT_CONFIG_JSONB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:

save/SaveV2/shape_and_slicesConst*A
value8B6B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:

save/SaveV2SaveV2
save/Constsave/SaveV2/tensor_namessave/SaveV2/shape_and_slicessave/Const_1save/Const_2save/Const_3dense/bias/Read/ReadVariableOp.training/Adam/dense/bias/m/Read/ReadVariableOp.training/Adam/dense/bias/v/Read/ReadVariableOp dense/kernel/Read/ReadVariableOp0training/Adam/dense/kernel/m/Read/ReadVariableOp0training/Adam/dense/kernel/v/Read/ReadVariableOpsave/Const_4 dense_1/bias/Read/ReadVariableOp0training/Adam/dense_1/bias/m/Read/ReadVariableOp0training/Adam/dense_1/bias/v/Read/ReadVariableOp"dense_1/kernel/Read/ReadVariableOp2training/Adam/dense_1/kernel/m/Read/ReadVariableOp2training/Adam/dense_1/kernel/v/Read/ReadVariableOpsave/Const_5beta_1/Read/ReadVariableOpbeta_2/Read/ReadVariableOpdecay/Read/ReadVariableOpepsilon/Read/ReadVariableOpiter/Read/ReadVariableOp!learning_rate/Read/ReadVariableOp*%
dtypes
2	
}
save/control_dependencyIdentity
save/Const^save/SaveV2*
T0*
_class
loc:@save/Const*
_output_shapes
: 
Ř
save/RestoreV2/tensor_namesConst"/device:CPU:0*ů

valueď
Bě
B/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-0/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-0/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB(optimizer/.ATTRIBUTES/OBJECT_CONFIG_JSONB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
Ł
save/RestoreV2/shape_and_slicesConst"/device:CPU:0*A
value8B6B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:

save/RestoreV2	RestoreV2
save/Constsave/RestoreV2/tensor_namessave/RestoreV2/shape_and_slices"/device:CPU:0*p
_output_shapes^
\:::::::::::::::::::::::*%
dtypes
2	

	save/NoOpNoOp

save/NoOp_1NoOp

save/NoOp_2NoOp
N
save/IdentityIdentitysave/RestoreV2:3*
_output_shapes
:*
T0
Q
save/AssignVariableOpAssignVariableOp
dense/biassave/Identity*
dtype0
P
save/Identity_1Identitysave/RestoreV2:4*
T0*
_output_shapes
:
e
save/AssignVariableOp_1AssignVariableOptraining/Adam/dense/bias/msave/Identity_1*
dtype0
P
save/Identity_2Identitysave/RestoreV2:5*
_output_shapes
:*
T0
e
save/AssignVariableOp_2AssignVariableOptraining/Adam/dense/bias/vsave/Identity_2*
dtype0
P
save/Identity_3Identitysave/RestoreV2:6*
_output_shapes
:*
T0
W
save/AssignVariableOp_3AssignVariableOpdense/kernelsave/Identity_3*
dtype0
P
save/Identity_4Identitysave/RestoreV2:7*
T0*
_output_shapes
:
g
save/AssignVariableOp_4AssignVariableOptraining/Adam/dense/kernel/msave/Identity_4*
dtype0
P
save/Identity_5Identitysave/RestoreV2:8*
T0*
_output_shapes
:
g
save/AssignVariableOp_5AssignVariableOptraining/Adam/dense/kernel/vsave/Identity_5*
dtype0

save/NoOp_3NoOp
Q
save/Identity_6Identitysave/RestoreV2:10*
T0*
_output_shapes
:
W
save/AssignVariableOp_6AssignVariableOpdense_1/biassave/Identity_6*
dtype0
Q
save/Identity_7Identitysave/RestoreV2:11*
_output_shapes
:*
T0
g
save/AssignVariableOp_7AssignVariableOptraining/Adam/dense_1/bias/msave/Identity_7*
dtype0
Q
save/Identity_8Identitysave/RestoreV2:12*
T0*
_output_shapes
:
g
save/AssignVariableOp_8AssignVariableOptraining/Adam/dense_1/bias/vsave/Identity_8*
dtype0
Q
save/Identity_9Identitysave/RestoreV2:13*
T0*
_output_shapes
:
Y
save/AssignVariableOp_9AssignVariableOpdense_1/kernelsave/Identity_9*
dtype0
R
save/Identity_10Identitysave/RestoreV2:14*
_output_shapes
:*
T0
k
save/AssignVariableOp_10AssignVariableOptraining/Adam/dense_1/kernel/msave/Identity_10*
dtype0
R
save/Identity_11Identitysave/RestoreV2:15*
_output_shapes
:*
T0
k
save/AssignVariableOp_11AssignVariableOptraining/Adam/dense_1/kernel/vsave/Identity_11*
dtype0

save/NoOp_4NoOp
R
save/Identity_12Identitysave/RestoreV2:17*
T0*
_output_shapes
:
S
save/AssignVariableOp_12AssignVariableOpbeta_1save/Identity_12*
dtype0
R
save/Identity_13Identitysave/RestoreV2:18*
T0*
_output_shapes
:
S
save/AssignVariableOp_13AssignVariableOpbeta_2save/Identity_13*
dtype0
R
save/Identity_14Identitysave/RestoreV2:19*
T0*
_output_shapes
:
R
save/AssignVariableOp_14AssignVariableOpdecaysave/Identity_14*
dtype0
R
save/Identity_15Identitysave/RestoreV2:20*
T0*
_output_shapes
:
T
save/AssignVariableOp_15AssignVariableOpepsilonsave/Identity_15*
dtype0
a
save/Identity_16Identitysave/RestoreV2:21"/device:CPU:0*
T0	*
_output_shapes
:
`
save/AssignVariableOp_16AssignVariableOpitersave/Identity_16"/device:CPU:0*
dtype0	
R
save/Identity_17Identitysave/RestoreV2:22*
T0*
_output_shapes
:
Z
save/AssignVariableOp_17AssignVariableOplearning_ratesave/Identity_17*
dtype0
 
save/restore_all/NoOpNoOp^save/AssignVariableOp^save/AssignVariableOp_1^save/AssignVariableOp_10^save/AssignVariableOp_11^save/AssignVariableOp_12^save/AssignVariableOp_13^save/AssignVariableOp_14^save/AssignVariableOp_15^save/AssignVariableOp_17^save/AssignVariableOp_2^save/AssignVariableOp_3^save/AssignVariableOp_4^save/AssignVariableOp_5^save/AssignVariableOp_6^save/AssignVariableOp_7^save/AssignVariableOp_8^save/AssignVariableOp_9
^save/NoOp^save/NoOp_1^save/NoOp_2^save/NoOp_3^save/NoOp_4
I
save/restore_all/NoOp_1NoOp^save/AssignVariableOp_16"/device:CPU:0
J
save/restore_allNoOp^save/restore_all/NoOp^save/restore_all/NoOp_1
0
init_1NoOp^count_1/Assign^total_1/Assign"D
save/Const:0save/control_dependency:0save/restore_all 5 @F8"ň
trainable_variablesÚ×
x
dense/kernel:0dense/kernel/Assign"dense/kernel/Read/ReadVariableOp:0(2)dense/kernel/Initializer/random_uniform:08
g
dense/bias:0dense/bias/Assign dense/bias/Read/ReadVariableOp:0(2dense/bias/Initializer/zeros:08

dense_1/kernel:0dense_1/kernel/Assign$dense_1/kernel/Read/ReadVariableOp:0(2+dense_1/kernel/Initializer/random_uniform:08
o
dense_1/bias:0dense_1/bias/Assign"dense_1/bias/Read/ReadVariableOp:0(2 dense_1/bias/Initializer/zeros:08"Í
local_variablesšś
Y
	total_1:0total_1/Assigntotal_1/Read/ReadVariableOp:0(2total_1/Initializer/zeros:0
Y
	count_1:0count_1/Assigncount_1/Read/ReadVariableOp:0(2count_1/Initializer/zeros:0" 
cond_contextřô
ć

rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/cond_textrloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_t:0 *
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:0
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1
rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_t:0č
rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0Ü
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1
Đr
tloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/cond_text_1rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_f:0*ú5
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Merge:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Merge:1
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:1
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:1
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:0
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:2
§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1
¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim:0
 loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims:0
Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
Ťloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1
Śloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim:0
˘loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1:0
Ąloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat:0
Śloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims:0
Ľloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const:0
Ľloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Shape:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0
rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_f:0
gloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank:0
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1:0č
rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
gloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch:0
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:02,
,
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/cond_textloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0 *ç(
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:0
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:2
§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1
¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim:0
 loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims:0
Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
Ťloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1
Śloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim:0
˘loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1:0
Ąloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat:0
Śloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims:0
Ľloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const:0
Ľloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Shape:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0Ô
§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0Ťloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1Ř
Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:12Ń
Î
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/cond_text_1loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0*

loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:1
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0¤
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0

oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/cond_textoloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t:0 *Ä
zloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency:0
oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0
ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t:0â
oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0

qloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/cond_text_1oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f:0*ź
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch:0
wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1:0
wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2:0
wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7:0
|loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1:0
oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0
ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f:0
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0
ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0ŕ
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3:0é
ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge:0uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch:0â
oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0ä
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1:0ă
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2:0"`
global_stepQO
M
iter:0iter/Assigniter/Read/ReadVariableOp:0(2iter/Initializer/zeros:0"Ž
	variables 
x
dense/kernel:0dense/kernel/Assign"dense/kernel/Read/ReadVariableOp:0(2)dense/kernel/Initializer/random_uniform:08
g
dense/bias:0dense/bias/Assign dense/bias/Read/ReadVariableOp:0(2dense/bias/Initializer/zeros:08

dense_1/kernel:0dense_1/kernel/Assign$dense_1/kernel/Read/ReadVariableOp:0(2+dense_1/kernel/Initializer/random_uniform:08
o
dense_1/bias:0dense_1/bias/Assign"dense_1/bias/Read/ReadVariableOp:0(2 dense_1/bias/Initializer/zeros:08
M
iter:0iter/Assigniter/Read/ReadVariableOp:0(2iter/Initializer/zeros:0
y
learning_rate:0learning_rate/Assign#learning_rate/Read/ReadVariableOp:0(2)learning_rate/Initializer/initial_value:0
Y
decay:0decay/Assigndecay/Read/ReadVariableOp:0(2!decay/Initializer/initial_value:0
]
beta_1:0beta_1/Assignbeta_1/Read/ReadVariableOp:0(2"beta_1/Initializer/initial_value:0
]
beta_2:0beta_2/Assignbeta_2/Read/ReadVariableOp:0(2"beta_2/Initializer/initial_value:0
a
	epsilon:0epsilon/Assignepsilon/Read/ReadVariableOp:0(2#epsilon/Initializer/initial_value:0
­
training/Adam/dense/kernel/m:0#training/Adam/dense/kernel/m/Assign2training/Adam/dense/kernel/m/Read/ReadVariableOp:0(20training/Adam/dense/kernel/m/Initializer/zeros:0
Ľ
training/Adam/dense/bias/m:0!training/Adam/dense/bias/m/Assign0training/Adam/dense/bias/m/Read/ReadVariableOp:0(2.training/Adam/dense/bias/m/Initializer/zeros:0
ľ
 training/Adam/dense_1/kernel/m:0%training/Adam/dense_1/kernel/m/Assign4training/Adam/dense_1/kernel/m/Read/ReadVariableOp:0(22training/Adam/dense_1/kernel/m/Initializer/zeros:0
­
training/Adam/dense_1/bias/m:0#training/Adam/dense_1/bias/m/Assign2training/Adam/dense_1/bias/m/Read/ReadVariableOp:0(20training/Adam/dense_1/bias/m/Initializer/zeros:0
­
training/Adam/dense/kernel/v:0#training/Adam/dense/kernel/v/Assign2training/Adam/dense/kernel/v/Read/ReadVariableOp:0(20training/Adam/dense/kernel/v/Initializer/zeros:0
Ľ
training/Adam/dense/bias/v:0!training/Adam/dense/bias/v/Assign0training/Adam/dense/bias/v/Read/ReadVariableOp:0(2.training/Adam/dense/bias/v/Initializer/zeros:0
ľ
 training/Adam/dense_1/kernel/v:0%training/Adam/dense_1/kernel/v/Assign4training/Adam/dense_1/kernel/v/Read/ReadVariableOp:0(22training/Adam/dense_1/kernel/v/Initializer/zeros:0
­
training/Adam/dense_1/bias/v:0#training/Adam/dense_1/bias/v/Assign2training/Adam/dense_1/bias/v/Read/ReadVariableOp:0(20training/Adam/dense_1/bias/v/Initializer/zeros:0*Q
__saved_model_train_op75
__saved_model_train_op
training_1/group_deps*ç
trainÝ
B
dense_1_target0
dense_1_target:0˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙
3
dense_input$
dense_input:0˙˙˙˙˙˙˙˙˙
9
metrics/accuracy/update_op
metric_op_wrapper:0 ,
metrics/accuracy/value
Identity_6:0 ?
predictions/dense_1(
dense_1/BiasAdd:0˙˙˙˙˙˙˙˙˙
loss

loss/mul:0 tensorflow/supervised/training*@
__saved_model_init_op'%
__saved_model_init_op
init_1¸°
Ó 
:
Add
x"T
y"T
z"T"
Ttype:
2	

ArgMax

input"T
	dimension"Tidx
output"output_type" 
Ttype:
2	"
Tidxtype0:
2	"
output_typetype0	:
2	
P
Assert
	condition
	
data2T"
T
list(type)(0"
	summarizeint
E
AssignAddVariableOp
resource
value"dtype"
dtypetype
B
AssignVariableOp
resource
value"dtype"
dtypetype
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
N
Cast	
x"SrcT	
y"DstT"
SrcTtype"
DstTtype"
Truncatebool( 
h
ConcatV2
values"T*N
axis"Tidx
output"T"
Nint(0"	
Ttype"
Tidxtype0:
2	
8
Const
output"dtype"
valuetensor"
dtypetype
š
DenseToDenseSetOperation	
set1"T	
set2"T
result_indices	
result_values"T
result_shape	"
set_operationstring"
validate_indicesbool("
Ttype:
	2	
5
DivNoNan
x"T
y"T
z"T"
Ttype:
2
B
Equal
x"T
y"T
z
"
Ttype:
2	

W

ExpandDims

input"T
dim"Tdim
output"T"	
Ttype"
Tdimtype0:
2	
^
Fill
dims"
index_type

value"T
output"T"	
Ttype"

index_typetype0:
2	
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	

Mean

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
N
Merge
inputs"T*N
output"T
value_index"	
Ttype"
Nint(0
=
Mul
x"T
y"T
z"T"
Ttype:
2	

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
X
PlaceholderWithDefault
input"dtype
output"dtype"
dtypetype"
shapeshape
~
RandomUniform

shape"T
output"dtype"
seedint "
seed2int "
dtypetype:
2"
Ttype:
2	
@
ReadVariableOp
resource
value"dtype"
dtypetype
E
Relu
features"T
activations"T"
Ttype:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
P
Shape

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
O
Size

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
G
SquaredDifference
x"T
y"T
z"T"
Ttype:

2	
:
Sub
x"T
y"T
z"T"
Ttype:
2	

Sum

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
M
Switch	
data"T
pred

output_false"T
output_true"T"	
Ttype
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape
9
VarIsInitializedOp
resource
is_initialized
"eval*2.0.0-alpha02v1.12.0-9492-g2c319fb4158ĺ
n
dense_inputPlaceholder*
shape:˙˙˙˙˙˙˙˙˙
*
dtype0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙


-dense/kernel/Initializer/random_uniform/shapeConst*
valueB"
      *
_class
loc:@dense/kernel*
dtype0*
_output_shapes
:

+dense/kernel/Initializer/random_uniform/minConst*
valueB
 *č!ż*
_class
loc:@dense/kernel*
dtype0*
_output_shapes
: 

+dense/kernel/Initializer/random_uniform/maxConst*
valueB
 *č!?*
_class
loc:@dense/kernel*
dtype0*
_output_shapes
: 
Ě
5dense/kernel/Initializer/random_uniform/RandomUniformRandomUniform-dense/kernel/Initializer/random_uniform/shape*
T0*
_class
loc:@dense/kernel*
dtype0*
_output_shapes

:

Î
+dense/kernel/Initializer/random_uniform/subSub+dense/kernel/Initializer/random_uniform/max+dense/kernel/Initializer/random_uniform/min*
_output_shapes
: *
T0*
_class
loc:@dense/kernel
ŕ
+dense/kernel/Initializer/random_uniform/mulMul5dense/kernel/Initializer/random_uniform/RandomUniform+dense/kernel/Initializer/random_uniform/sub*
T0*
_class
loc:@dense/kernel*
_output_shapes

:

Ň
'dense/kernel/Initializer/random_uniformAdd+dense/kernel/Initializer/random_uniform/mul+dense/kernel/Initializer/random_uniform/min*
T0*
_class
loc:@dense/kernel*
_output_shapes

:


dense/kernelVarHandleOp*
shape
:
*
shared_namedense/kernel*
_class
loc:@dense/kernel*
dtype0*
_output_shapes
: 
i
-dense/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense/kernel*
_output_shapes
: 

dense/kernel/AssignAssignVariableOpdense/kernel'dense/kernel/Initializer/random_uniform*
_class
loc:@dense/kernel*
dtype0

 dense/kernel/Read/ReadVariableOpReadVariableOpdense/kernel*
_class
loc:@dense/kernel*
dtype0*
_output_shapes

:


dense/bias/Initializer/zerosConst*
valueB*    *
_class
loc:@dense/bias*
dtype0*
_output_shapes
:


dense/biasVarHandleOp*
shape:*
shared_name
dense/bias*
_class
loc:@dense/bias*
dtype0*
_output_shapes
: 
e
+dense/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOp
dense/bias*
_output_shapes
: 
{
dense/bias/AssignAssignVariableOp
dense/biasdense/bias/Initializer/zeros*
_class
loc:@dense/bias*
dtype0

dense/bias/Read/ReadVariableOpReadVariableOp
dense/bias*
_class
loc:@dense/bias*
dtype0*
_output_shapes
:
h
dense/MatMul/ReadVariableOpReadVariableOpdense/kernel*
dtype0*
_output_shapes

:

r
dense/MatMulMatMuldense_inputdense/MatMul/ReadVariableOp*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
c
dense/BiasAdd/ReadVariableOpReadVariableOp
dense/bias*
dtype0*
_output_shapes
:
v
dense/BiasAddBiasAdddense/MatMuldense/BiasAdd/ReadVariableOp*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
S

dense/ReluReludense/BiasAdd*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
Ł
/dense_1/kernel/Initializer/random_uniform/shapeConst*
dtype0*
_output_shapes
:*
valueB"      *!
_class
loc:@dense_1/kernel

-dense_1/kernel/Initializer/random_uniform/minConst*
valueB
 *mż*!
_class
loc:@dense_1/kernel*
dtype0*
_output_shapes
: 

-dense_1/kernel/Initializer/random_uniform/maxConst*
valueB
 *m?*!
_class
loc:@dense_1/kernel*
dtype0*
_output_shapes
: 
Ň
7dense_1/kernel/Initializer/random_uniform/RandomUniformRandomUniform/dense_1/kernel/Initializer/random_uniform/shape*
T0*!
_class
loc:@dense_1/kernel*
dtype0*
_output_shapes

:
Ö
-dense_1/kernel/Initializer/random_uniform/subSub-dense_1/kernel/Initializer/random_uniform/max-dense_1/kernel/Initializer/random_uniform/min*
T0*!
_class
loc:@dense_1/kernel*
_output_shapes
: 
č
-dense_1/kernel/Initializer/random_uniform/mulMul7dense_1/kernel/Initializer/random_uniform/RandomUniform-dense_1/kernel/Initializer/random_uniform/sub*
T0*!
_class
loc:@dense_1/kernel*
_output_shapes

:
Ú
)dense_1/kernel/Initializer/random_uniformAdd-dense_1/kernel/Initializer/random_uniform/mul-dense_1/kernel/Initializer/random_uniform/min*
T0*!
_class
loc:@dense_1/kernel*
_output_shapes

:

dense_1/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape
:*
shared_namedense_1/kernel*!
_class
loc:@dense_1/kernel
m
/dense_1/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_1/kernel*
_output_shapes
: 

dense_1/kernel/AssignAssignVariableOpdense_1/kernel)dense_1/kernel/Initializer/random_uniform*!
_class
loc:@dense_1/kernel*
dtype0

"dense_1/kernel/Read/ReadVariableOpReadVariableOpdense_1/kernel*!
_class
loc:@dense_1/kernel*
dtype0*
_output_shapes

:

dense_1/bias/Initializer/zerosConst*
valueB*    *
_class
loc:@dense_1/bias*
dtype0*
_output_shapes
:

dense_1/biasVarHandleOp*
dtype0*
_output_shapes
: *
shape:*
shared_namedense_1/bias*
_class
loc:@dense_1/bias
i
-dense_1/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_1/bias*
_output_shapes
: 

dense_1/bias/AssignAssignVariableOpdense_1/biasdense_1/bias/Initializer/zeros*
_class
loc:@dense_1/bias*
dtype0

 dense_1/bias/Read/ReadVariableOpReadVariableOpdense_1/bias*
dtype0*
_output_shapes
:*
_class
loc:@dense_1/bias
l
dense_1/MatMul/ReadVariableOpReadVariableOpdense_1/kernel*
dtype0*
_output_shapes

:
u
dense_1/MatMulMatMul
dense/Reludense_1/MatMul/ReadVariableOp*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
g
dense_1/BiasAdd/ReadVariableOpReadVariableOpdense_1/bias*
dtype0*
_output_shapes
:
|
dense_1/BiasAddBiasAdddense_1/MatMuldense_1/BiasAdd/ReadVariableOp*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0

dense_1_targetPlaceholder*
dtype0*0
_output_shapes
:˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙*%
shape:˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙
R
ConstConst*
dtype0*
_output_shapes
:*
valueB*  ?

dense_1_sample_weightsPlaceholderWithDefaultConst*
dtype0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*
shape:˙˙˙˙˙˙˙˙˙
v
total/Initializer/zerosConst*
valueB
 *    *
_class

loc:@total*
dtype0*
_output_shapes
: 
x
totalVarHandleOp*
shared_nametotal*
_class

loc:@total*
dtype0*
_output_shapes
: *
shape: 
[
&total/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal*
_output_shapes
: 
g
total/AssignAssignVariableOptotaltotal/Initializer/zeros*
_class

loc:@total*
dtype0
q
total/Read/ReadVariableOpReadVariableOptotal*
_class

loc:@total*
dtype0*
_output_shapes
: 
v
count/Initializer/zerosConst*
valueB
 *    *
_class

loc:@count*
dtype0*
_output_shapes
: 
x
countVarHandleOp*
_class

loc:@count*
dtype0*
_output_shapes
: *
shape: *
shared_namecount
[
&count/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount*
_output_shapes
: 
g
count/AssignAssignVariableOpcountcount/Initializer/zeros*
_class

loc:@count*
dtype0
q
count/Read/ReadVariableOpReadVariableOpcount*
dtype0*
_output_shapes
: *
_class

loc:@count
l
!metrics/accuracy/ArgMax/dimensionConst*
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0*
_output_shapes
: 

metrics/accuracy/ArgMaxArgMaxdense_1_target!metrics/accuracy/ArgMax/dimension*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
n
#metrics/accuracy/ArgMax_1/dimensionConst*
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0*
_output_shapes
: 

metrics/accuracy/ArgMax_1ArgMaxdense_1/BiasAdd#metrics/accuracy/ArgMax_1/dimension*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙

metrics/accuracy/EqualEqualmetrics/accuracy/ArgMaxmetrics/accuracy/ArgMax_1*
T0	*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
r
metrics/accuracy/CastCastmetrics/accuracy/Equal*

SrcT0
*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*

DstT0
`
metrics/accuracy/ConstConst*
valueB: *
dtype0*
_output_shapes
:
k
metrics/accuracy/SumSummetrics/accuracy/Castmetrics/accuracy/Const*
_output_shapes
: *
T0
e
$metrics/accuracy/AssignAddVariableOpAssignAddVariableOptotalmetrics/accuracy/Sum*
dtype0

metrics/accuracy/ReadVariableOpReadVariableOptotal%^metrics/accuracy/AssignAddVariableOp^metrics/accuracy/Sum*
dtype0*
_output_shapes
: 
U
metrics/accuracy/SizeSizemetrics/accuracy/Cast*
_output_shapes
: *
T0
f
metrics/accuracy/Cast_1Castmetrics/accuracy/Size*
_output_shapes
: *

DstT0*

SrcT0

&metrics/accuracy/AssignAddVariableOp_1AssignAddVariableOpcountmetrics/accuracy/Cast_1%^metrics/accuracy/AssignAddVariableOp*
dtype0
Ż
!metrics/accuracy/ReadVariableOp_1ReadVariableOpcount%^metrics/accuracy/AssignAddVariableOp'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

*metrics/accuracy/div_no_nan/ReadVariableOpReadVariableOptotal'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 

,metrics/accuracy/div_no_nan/ReadVariableOp_1ReadVariableOpcount'^metrics/accuracy/AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
˘
metrics/accuracy/div_no_nanDivNoNan*metrics/accuracy/div_no_nan/ReadVariableOp,metrics/accuracy/div_no_nan/ReadVariableOp_1*
_output_shapes
: *
T0
c
metrics/accuracy/IdentityIdentitymetrics/accuracy/div_no_nan*
T0*
_output_shapes
: 
n
#metrics/accuracy/ArgMax_2/dimensionConst*
_output_shapes
: *
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0

metrics/accuracy/ArgMax_2ArgMaxdense_1_target#metrics/accuracy/ArgMax_2/dimension*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
n
#metrics/accuracy/ArgMax_3/dimensionConst*
_output_shapes
: *
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0

metrics/accuracy/ArgMax_3ArgMaxdense_1/BiasAdd#metrics/accuracy/ArgMax_3/dimension*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0

metrics/accuracy/Equal_1Equalmetrics/accuracy/ArgMax_2metrics/accuracy/ArgMax_3*
T0	*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
v
metrics/accuracy/Cast_2Castmetrics/accuracy/Equal_1*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*

DstT0*

SrcT0

b
metrics/accuracy/Const_1Const*
dtype0*
_output_shapes
:*
valueB: 
q
metrics/accuracy/MeanMeanmetrics/accuracy/Cast_2metrics/accuracy/Const_1*
T0*
_output_shapes
: 

6loss/dense_1_loss/mean_squared_error/SquaredDifferenceSquaredDifferencedense_1/BiasAdddense_1_target*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙

;loss/dense_1_loss/mean_squared_error/Mean/reduction_indicesConst*
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0*
_output_shapes
: 
Ô
)loss/dense_1_loss/mean_squared_error/MeanMean6loss/dense_1_loss/mean_squared_error/SquaredDifference;loss/dense_1_loss/mean_squared_error/Mean/reduction_indices*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
­
gloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapeShapedense_1_sample_weights*
T0*
_output_shapes
:
¨
floss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rankConst*
value	B :*
dtype0*
_output_shapes
: 
ż
floss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapeShape)loss/dense_1_loss/mean_squared_error/Mean*
_output_shapes
:*
T0
§
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rankConst*
value	B :*
dtype0*
_output_shapes
: 
§
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar/xConst*
value	B : *
dtype0*
_output_shapes
: 
Ü
closs/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarEqualeloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar/xfloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank*
T0*
_output_shapes
: 
ć
oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/SwitchSwitchcloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarcloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar*
T0
*
_output_shapes
: : 

qloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_tIdentityqloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch:1*
T0
*
_output_shapes
: 

qloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_fIdentityoloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch*
T0
*
_output_shapes
: 

ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_idIdentitycloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar*
T0
*
_output_shapes
: 
í
qloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1Switchcloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
_output_shapes
: : *
T0
*v
_classl
jhloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar
î
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankEqualloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switchloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1*
T0*
_output_shapes
: 

loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/SwitchSwitcheloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rankploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
T0*x
_classn
ljloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank*
_output_shapes
: : 

loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1Switchfloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rankploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*
T0*y
_classo
mkloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank*
_output_shapes
: : 
Ű
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/SwitchSwitchloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
T0
*
_output_shapes
: : 
Ç
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_tIdentityloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:1*
_output_shapes
: *
T0

Ĺ
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_fIdentityloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch*
T0
*
_output_shapes
: 
Ę
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_idIdentityloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank*
T0
*
_output_shapes
: 
ý
˘loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dimConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0*
_output_shapes
: 
§
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims
ExpandDimsŠloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1˘loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim*
T0*
_output_shapes

:
°
Ľloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/SwitchSwitchfloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapeploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id*y
_classo
mkloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::*
T0

§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1SwitchĽloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switchloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0*y
_classo
mkloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::

Łloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/ShapeConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
_output_shapes
:*
valueB"      *
dtype0
ő
Łloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/ConstConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
value	B :*
dtype0*
_output_shapes
: 

loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_likeFillŁloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/ShapeŁloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const*
T0*
_output_shapes

:
ń
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axisConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
value	B :*
dtype0*
_output_shapes
: 
ź
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concatConcatV2loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDimsloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_likeloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis*
T0*
N*
_output_shapes

:
˙
¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dimConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
valueB :
˙˙˙˙˙˙˙˙˙*
dtype0*
_output_shapes
: 
­
 loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1
ExpandDimsŤloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim*
_output_shapes

:*
T0
´
§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/SwitchSwitchgloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapeploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id* 
_output_shapes
::*
T0*z
_classp
nlloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape

Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1Switch§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switchloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
T0*z
_classp
nlloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::
č
Źloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperationDenseToDenseSetOperation loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat*
T0*<
_output_shapes*
(:˙˙˙˙˙˙˙˙˙:˙˙˙˙˙˙˙˙˙:*
set_operationa-b
˙
¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dimsSizeŽloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1*
_output_shapes
: *
T0
ç
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/xConst^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t*
_output_shapes
: *
value	B : *
dtype0
ý
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dimsEqualloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims*
T0*
_output_shapes
: 

loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1Switchloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rankloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id*
_output_shapes
: : *
T0
*Ľ
_class
loc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank
â
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/MergeMergeloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims*
N*
_output_shapes
: : *
T0

Ł
nloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/MergeMergeloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Mergesloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1*
T0
*
N*
_output_shapes
: : 
Ç
_loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/ConstConst*8
value/B- B'weights can not be broadcast to values.*
dtype0*
_output_shapes
: 
°
aloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_1Const*
valueB Bweights.shape=*
dtype0*
_output_shapes
: 
ş
aloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_2Const*)
value B Bdense_1_sample_weights:0*
dtype0*
_output_shapes
: 
Ż
aloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_3Const*
dtype0*
_output_shapes
: *
valueB Bvalues.shape=
Í
aloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_4Const*<
value3B1 B+loss/dense_1_loss/mean_squared_error/Mean:0*
dtype0*
_output_shapes
: 
Ź
aloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/Const_5Const*
valueB B
is_scalar=*
dtype0*
_output_shapes
: 
ů
lloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/SwitchSwitchnloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Mergenloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge*
T0
*
_output_shapes
: : 

nloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_tIdentitynloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Switch:1*
_output_shapes
: *
T0


nloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_fIdentitylloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Switch*
T0
*
_output_shapes
: 

mloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_idIdentitynloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge*
_output_shapes
: *
T0

ă
jloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/NoOpNoOpo^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t

xloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependencyIdentitynloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_tk^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/NoOp*
T0
*
_classw
usloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t*
_output_shapes
: 
Ě
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
dtype0*
_output_shapes
: *8
value/B- B'weights can not be broadcast to values.
ł
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
dtype0*
_output_shapes
: *
valueB Bweights.shape=
˝
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*)
value B Bdense_1_sample_weights:0*
dtype0*
_output_shapes
: 
˛
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
valueB Bvalues.shape=*
dtype0*
_output_shapes
: 
Đ
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*<
value3B1 B+loss/dense_1_loss/mean_squared_error/Mean:0*
dtype0*
_output_shapes
: 
Ż
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7Consto^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f*
valueB B
is_scalar=*
dtype0*
_output_shapes
: 
˘

lloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/AssertAssertsloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switchsloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3*
T
2	


sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/SwitchSwitchnloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Mergemloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0
*
_classw
usloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge*
_output_shapes
: : 
ţ
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1Switchgloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shapemloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0*z
_classp
nlloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape* 
_output_shapes
::
ü
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2Switchfloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shapemloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
T0*y
_classo
mkloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape* 
_output_shapes
::
î
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3Switchcloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalarmloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id*
_output_shapes
: : *
T0
*v
_classl
jhloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar

zloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1Identitynloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_fm^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert*
_output_shapes
: *
T0
*
_classw
usloc:@loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f

kloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/MergeMergezloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1xloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency*
T0
*
N*
_output_shapes
: : 

Tloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ShapeShape)loss/dense_1_loss/mean_squared_error/Meanl^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Merge*
T0*
_output_shapes
:

Tloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ConstConstl^loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Merge*
dtype0*
_output_shapes
: *
valueB
 *  ?
°
Nloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_likeFillTloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/ShapeTloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like/Const*#
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
á
Dloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weightsMuldense_1_sample_weightsNloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/ones_like*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
Ü
6loss/dense_1_loss/mean_squared_error/weighted_loss/MulMul)loss/dense_1_loss/mean_squared_error/MeanDloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights*
T0*#
_output_shapes
:˙˙˙˙˙˙˙˙˙
a
loss/dense_1_loss/ConstConst*
valueB: *
dtype0*
_output_shapes
:

loss/dense_1_loss/SumSum6loss/dense_1_loss/mean_squared_error/weighted_loss/Mulloss/dense_1_loss/Const*
_output_shapes
: *
T0

loss/dense_1_loss/num_elementsSize6loss/dense_1_loss/mean_squared_error/weighted_loss/Mul*
T0*
_output_shapes
: 
{
#loss/dense_1_loss/num_elements/CastCastloss/dense_1_loss/num_elements*

SrcT0*
_output_shapes
: *

DstT0
\
loss/dense_1_loss/mul/xConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 
{
loss/dense_1_loss/mulMulloss/dense_1_loss/mul/x#loss/dense_1_loss/num_elements/Cast*
_output_shapes
: *
T0
\
loss/dense_1_loss/Const_1Const*
valueB *
dtype0*
_output_shapes
: 
q
loss/dense_1_loss/Sum_1Sumloss/dense_1_loss/Sumloss/dense_1_loss/Const_1*
_output_shapes
: *
T0
t
loss/dense_1_loss/valueDivNoNanloss/dense_1_loss/Sum_1loss/dense_1_loss/mul*
T0*
_output_shapes
: 
O

loss/mul/xConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 
U
loss/mulMul
loss/mul/xloss/dense_1_loss/value*
_output_shapes
: *
T0
q
iter/Initializer/zerosConst*
_output_shapes
: *
value	B	 R *
_class
	loc:@iter*
dtype0	

iterVarHandleOp"/device:CPU:0*
shared_nameiter*
_class
	loc:@iter*
dtype0	*
_output_shapes
: *
shape: 
h
%iter/IsInitialized/VarIsInitializedOpVarIsInitializedOpiter"/device:CPU:0*
_output_shapes
: 
r
iter/AssignAssignVariableOpiteriter/Initializer/zeros"/device:CPU:0*
_class
	loc:@iter*
dtype0	
}
iter/Read/ReadVariableOpReadVariableOpiter"/device:CPU:0*
_class
	loc:@iter*
dtype0	*
_output_shapes
: 

'learning_rate/Initializer/initial_valueConst*
valueB
 *o:* 
_class
loc:@learning_rate*
dtype0*
_output_shapes
: 

learning_rateVarHandleOp*
shared_namelearning_rate* 
_class
loc:@learning_rate*
dtype0*
_output_shapes
: *
shape: 
k
.learning_rate/IsInitialized/VarIsInitializedOpVarIsInitializedOplearning_rate*
_output_shapes
: 

learning_rate/AssignAssignVariableOplearning_rate'learning_rate/Initializer/initial_value* 
_class
loc:@learning_rate*
dtype0

!learning_rate/Read/ReadVariableOpReadVariableOplearning_rate* 
_class
loc:@learning_rate*
dtype0*
_output_shapes
: 
~
decay/Initializer/initial_valueConst*
_output_shapes
: *
valueB
 *    *
_class

loc:@decay*
dtype0
x
decayVarHandleOp*
shared_namedecay*
_class

loc:@decay*
dtype0*
_output_shapes
: *
shape: 
[
&decay/IsInitialized/VarIsInitializedOpVarIsInitializedOpdecay*
_output_shapes
: 
o
decay/AssignAssignVariableOpdecaydecay/Initializer/initial_value*
_class

loc:@decay*
dtype0
q
decay/Read/ReadVariableOpReadVariableOpdecay*
_class

loc:@decay*
dtype0*
_output_shapes
: 

 beta_1/Initializer/initial_valueConst*
_output_shapes
: *
valueB
 *fff?*
_class
loc:@beta_1*
dtype0
{
beta_1VarHandleOp*
shape: *
shared_namebeta_1*
_class
loc:@beta_1*
dtype0*
_output_shapes
: 
]
'beta_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpbeta_1*
_output_shapes
: 
s
beta_1/AssignAssignVariableOpbeta_1 beta_1/Initializer/initial_value*
_class
loc:@beta_1*
dtype0
t
beta_1/Read/ReadVariableOpReadVariableOpbeta_1*
dtype0*
_output_shapes
: *
_class
loc:@beta_1

 beta_2/Initializer/initial_valueConst*
_output_shapes
: *
valueB
 *wž?*
_class
loc:@beta_2*
dtype0
{
beta_2VarHandleOp*
shape: *
shared_namebeta_2*
_class
loc:@beta_2*
dtype0*
_output_shapes
: 
]
'beta_2/IsInitialized/VarIsInitializedOpVarIsInitializedOpbeta_2*
_output_shapes
: 
s
beta_2/AssignAssignVariableOpbeta_2 beta_2/Initializer/initial_value*
_class
loc:@beta_2*
dtype0
t
beta_2/Read/ReadVariableOpReadVariableOpbeta_2*
_class
loc:@beta_2*
dtype0*
_output_shapes
: 

!epsilon/Initializer/initial_valueConst*
valueB
 *żÖ3*
_class
loc:@epsilon*
dtype0*
_output_shapes
: 
~
epsilonVarHandleOp*
_class
loc:@epsilon*
dtype0*
_output_shapes
: *
shape: *
shared_name	epsilon
_
(epsilon/IsInitialized/VarIsInitializedOpVarIsInitializedOpepsilon*
_output_shapes
: 
w
epsilon/AssignAssignVariableOpepsilon!epsilon/Initializer/initial_value*
_class
loc:@epsilon*
dtype0
w
epsilon/Read/ReadVariableOpReadVariableOpepsilon*
dtype0*
_output_shapes
: *
_class
loc:@epsilon
@
evaluation/group_depsNoOp	^loss/mul^metrics/accuracy/Mean
W
Const_1Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_2Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_3Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_4Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
G
VarIsInitializedOpVarIsInitializedOptotal*
_output_shapes
: 
H
VarIsInitializedOp_1VarIsInitializedOpiter*
_output_shapes
: 
J
VarIsInitializedOp_2VarIsInitializedOpbeta_1*
_output_shapes
: 
Q
VarIsInitializedOp_3VarIsInitializedOplearning_rate*
_output_shapes
: 
I
VarIsInitializedOp_4VarIsInitializedOpcount*
_output_shapes
: 
J
VarIsInitializedOp_5VarIsInitializedOpbeta_2*
_output_shapes
: 
K
VarIsInitializedOp_6VarIsInitializedOpepsilon*
_output_shapes
: 
P
VarIsInitializedOp_7VarIsInitializedOpdense/kernel*
_output_shapes
: 
R
VarIsInitializedOp_8VarIsInitializedOpdense_1/kernel*
_output_shapes
: 
N
VarIsInitializedOp_9VarIsInitializedOp
dense/bias*
_output_shapes
: 
Q
VarIsInitializedOp_10VarIsInitializedOpdense_1/bias*
_output_shapes
: 
J
VarIsInitializedOp_11VarIsInitializedOpdecay*
_output_shapes
: 
Ţ
	init/NoOpNoOp^beta_1/Assign^beta_2/Assign^count/Assign^decay/Assign^dense/bias/Assign^dense/kernel/Assign^dense_1/bias/Assign^dense_1/kernel/Assign^epsilon/Assign^learning_rate/Assign^total/Assign
0
init/NoOp_1NoOp^iter/Assign"/device:CPU:0
&
initNoOp
^init/NoOp^init/NoOp_1
W
Const_5Const"/device:CPU:0*
_output_shapes
: *
valueB B *
dtype0
\
Const_6Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB Bmodel
W
Const_7Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_8Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_9Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
X
Const_10Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
X
Const_11Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
¤
RestoreV2/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
r
RestoreV2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

	RestoreV2	RestoreV2Const_6RestoreV2/tensor_namesRestoreV2/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
B
IdentityIdentity	RestoreV2*
_output_shapes
:*
T0
I
AssignVariableOpAssignVariableOpdense/kernelIdentity*
dtype0
¤
RestoreV2_1/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_1	RestoreV2Const_6RestoreV2_1/tensor_namesRestoreV2_1/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_1IdentityRestoreV2_1*
T0*
_output_shapes
:
K
AssignVariableOp_1AssignVariableOp
dense/bias
Identity_1*
dtype0
Ś
RestoreV2_2/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_2	RestoreV2Const_6RestoreV2_2/tensor_namesRestoreV2_2/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_2IdentityRestoreV2_2*
_output_shapes
:*
T0
O
AssignVariableOp_2AssignVariableOpdense_1/kernel
Identity_2*
dtype0
¤
RestoreV2_3/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_3/shape_and_slicesConst"/device:CPU:0*
dtype0*
_output_shapes
:*
valueB
B 

RestoreV2_3	RestoreV2Const_6RestoreV2_3/tensor_namesRestoreV2_3/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_3IdentityRestoreV2_3*
_output_shapes
:*
T0
M
AssignVariableOp_3AssignVariableOpdense_1/bias
Identity_3*
dtype0

RestoreV2_4/tensor_namesConst"/device:CPU:0*>
value5B3B)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_4/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_4	RestoreV2Const_6RestoreV2_4/tensor_namesRestoreV2_4/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2	
U

Identity_4IdentityRestoreV2_4"/device:CPU:0*
_output_shapes
:*
T0	
T
AssignVariableOp_4AssignVariableOpiter
Identity_4"/device:CPU:0*
dtype0	
˘
RestoreV2_5/tensor_namesConst"/device:CPU:0*G
value>B<B2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_5/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_5	RestoreV2Const_6RestoreV2_5/tensor_namesRestoreV2_5/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_5IdentityRestoreV2_5*
_output_shapes
:*
T0
N
AssignVariableOp_5AssignVariableOplearning_rate
Identity_5*
dtype0

RestoreV2_6/tensor_namesConst"/device:CPU:0*?
value6B4B*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_6/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
valueB
B *
dtype0

RestoreV2_6	RestoreV2Const_6RestoreV2_6/tensor_namesRestoreV2_6/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_6IdentityRestoreV2_6*
T0*
_output_shapes
:
F
AssignVariableOp_6AssignVariableOpdecay
Identity_6*
dtype0

RestoreV2_7/tensor_namesConst"/device:CPU:0*
_output_shapes
:*@
value7B5B+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE*
dtype0
t
RestoreV2_7/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_7	RestoreV2Const_6RestoreV2_7/tensor_namesRestoreV2_7/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_7IdentityRestoreV2_7*
T0*
_output_shapes
:
G
AssignVariableOp_7AssignVariableOpbeta_1
Identity_7*
dtype0

RestoreV2_8/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*@
value7B5B+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE
t
RestoreV2_8/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_8	RestoreV2Const_6RestoreV2_8/tensor_namesRestoreV2_8/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_8IdentityRestoreV2_8*
T0*
_output_shapes
:
G
AssignVariableOp_8AssignVariableOpbeta_2
Identity_8*
dtype0

RestoreV2_9/tensor_namesConst"/device:CPU:0*A
value8B6B,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_9/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_9	RestoreV2Const_6RestoreV2_9/tensor_namesRestoreV2_9/shape_and_slices"/device:CPU:0*
dtypes
2*
_output_shapes
:
F

Identity_9IdentityRestoreV2_9*
T0*
_output_shapes
:
H
AssignVariableOp_9AssignVariableOpepsilon
Identity_9*
dtype0
z
total_1/Initializer/zerosConst*
valueB
 *    *
_class
loc:@total_1*
dtype0*
_output_shapes
: 
~
total_1VarHandleOp*
shape: *
shared_name	total_1*
_class
loc:@total_1*
dtype0*
_output_shapes
: 
_
(total_1/IsInitialized/VarIsInitializedOpVarIsInitializedOptotal_1*
_output_shapes
: 
o
total_1/AssignAssignVariableOptotal_1total_1/Initializer/zeros*
_class
loc:@total_1*
dtype0
w
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_class
loc:@total_1*
dtype0*
_output_shapes
: 
z
count_1/Initializer/zerosConst*
valueB
 *    *
_class
loc:@count_1*
dtype0*
_output_shapes
: 
~
count_1VarHandleOp*
shared_name	count_1*
_class
loc:@count_1*
dtype0*
_output_shapes
: *
shape: 
_
(count_1/IsInitialized/VarIsInitializedOpVarIsInitializedOpcount_1*
_output_shapes
: 
o
count_1/AssignAssignVariableOpcount_1count_1/Initializer/zeros*
_class
loc:@count_1*
dtype0
w
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_class
loc:@count_1*
dtype0*
_output_shapes
: 
K
Const_12Const*
valueB *
dtype0*
_output_shapes
: 
L
SumSummetrics/accuracy/MeanConst_12*
_output_shapes
: *
T0
E
AssignAddVariableOpAssignAddVariableOptotal_1Sum*
dtype0
j
ReadVariableOpReadVariableOptotal_1^AssignAddVariableOp^Sum*
dtype0*
_output_shapes
: 
F
SizeConst*
_output_shapes
: *
value	B :*
dtype0
B
CastCastSize*
_output_shapes
: *

DstT0*

SrcT0
^
AssignAddVariableOp_1AssignAddVariableOpcount_1Cast^AssignAddVariableOp*
dtype0
~
ReadVariableOp_1ReadVariableOpcount_1^AssignAddVariableOp^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
q
div_no_nan/ReadVariableOpReadVariableOptotal_1^AssignAddVariableOp_1*
_output_shapes
: *
dtype0
s
div_no_nan/ReadVariableOp_1ReadVariableOpcount_1^AssignAddVariableOp_1*
dtype0*
_output_shapes
: 
o

div_no_nanDivNoNandiv_no_nan/ReadVariableOpdiv_no_nan/ReadVariableOp_1*
_output_shapes
: *
T0
D
Identity_10Identity
div_no_nan*
T0*
_output_shapes
: 
[
div_no_nan_1/ReadVariableOpReadVariableOptotal_1*
dtype0*
_output_shapes
: 
]
div_no_nan_1/ReadVariableOp_1ReadVariableOpcount_1*
dtype0*
_output_shapes
: 
u
div_no_nan_1DivNoNandiv_no_nan_1/ReadVariableOpdiv_no_nan_1/ReadVariableOp_1*
T0*
_output_shapes
: 
F
Identity_11Identitydiv_no_nan_1*
T0*
_output_shapes
: 
l
metric_op_wrapperConst^AssignAddVariableOp_1*
valueB *
dtype0*
_output_shapes
: 
Y
save/filename/inputConst*
valueB Bmodel*
dtype0*
_output_shapes
: 
n
save/filenamePlaceholderWithDefaultsave/filename/input*
dtype0*
_output_shapes
: *
shape: 
e

save/ConstPlaceholderWithDefaultsave/filename*
dtype0*
_output_shapes
: *
shape: 

save/Const_1Const*Í
valueĂBŔ Bš{"class_name": "Sequential", "config": {"layers": [{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "batch_input_shape": [null, 10], "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense", "trainable": true, "units": 5, "use_bias": true}}, {"class_name": "Dense", "config": {"activation": "linear", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense_1", "trainable": true, "units": 2, "use_bias": true}}], "name": "sequential"}}*
dtype0*
_output_shapes
: 
Ö
save/Const_2Const*
valueB B{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 10], "dtype": "float32", "name": "dense_input", "sparse": false}}*
dtype0*
_output_shapes
: 

save/Const_3Const*Ř
valueÎBË BÄ{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "batch_input_shape": [null, 10], "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense", "trainable": true, "units": 5, "use_bias": true}}*
dtype0*
_output_shapes
: 
ř
save/Const_4Const*
dtype0*
_output_shapes
: *ť
valueąBŽ B§{"class_name": "Dense", "config": {"activation": "linear", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense_1", "trainable": true, "units": 2, "use_bias": true}}
N
save/VarIsInitializedOpVarIsInitializedOptotal_1*
_output_shapes
: 
P
save/VarIsInitializedOp_1VarIsInitializedOpcount_1*
_output_shapes
: 
3
	save/initNoOp^count_1/Assign^total_1/Assign
Ş
save/Const_5Const*
dtype0*
_output_shapes
: *í
valueăBŕ BŮ{"class_name": "Adam", "config": {"amsgrad": false, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "decay": 0.0, "epsilon": 1.0000000116860974e-07, "learning_rate": 0.0010000000474974513, "name": "Adam"}}
Ž
save/SaveV2/tensor_namesConst*á
value×BÔB/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB(optimizer/.ATTRIBUTES/OBJECT_CONFIG_JSONB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:

save/SaveV2/shape_and_slicesConst*1
value(B&B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:
ó
save/SaveV2SaveV2
save/Constsave/SaveV2/tensor_namessave/SaveV2/shape_and_slicessave/Const_1save/Const_2save/Const_3dense/bias/Read/ReadVariableOp dense/kernel/Read/ReadVariableOpsave/Const_4 dense_1/bias/Read/ReadVariableOp"dense_1/kernel/Read/ReadVariableOpsave/Const_5beta_1/Read/ReadVariableOpbeta_2/Read/ReadVariableOpdecay/Read/ReadVariableOpepsilon/Read/ReadVariableOpiter/Read/ReadVariableOp!learning_rate/Read/ReadVariableOp*
dtypes
2	
}
save/control_dependencyIdentity
save/Const^save/SaveV2*
T0*
_class
loc:@save/Const*
_output_shapes
: 
Ŕ
save/RestoreV2/tensor_namesConst"/device:CPU:0*
dtype0*
_output_shapes
:*á
value×BÔB/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB(optimizer/.ATTRIBUTES/OBJECT_CONFIG_JSONB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB,optimizer/epsilon/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE

save/RestoreV2/shape_and_slicesConst"/device:CPU:0*1
value(B&B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:
ĺ
save/RestoreV2	RestoreV2
save/Constsave/RestoreV2/tensor_namessave/RestoreV2/shape_and_slices"/device:CPU:0*P
_output_shapes>
<:::::::::::::::*
dtypes
2	

	save/NoOpNoOp

save/NoOp_1NoOp

save/NoOp_2NoOp
N
save/IdentityIdentitysave/RestoreV2:3*
T0*
_output_shapes
:
Q
save/AssignVariableOpAssignVariableOp
dense/biassave/Identity*
dtype0
P
save/Identity_1Identitysave/RestoreV2:4*
T0*
_output_shapes
:
W
save/AssignVariableOp_1AssignVariableOpdense/kernelsave/Identity_1*
dtype0

save/NoOp_3NoOp
P
save/Identity_2Identitysave/RestoreV2:6*
T0*
_output_shapes
:
W
save/AssignVariableOp_2AssignVariableOpdense_1/biassave/Identity_2*
dtype0
P
save/Identity_3Identitysave/RestoreV2:7*
_output_shapes
:*
T0
Y
save/AssignVariableOp_3AssignVariableOpdense_1/kernelsave/Identity_3*
dtype0

save/NoOp_4NoOp
P
save/Identity_4Identitysave/RestoreV2:9*
T0*
_output_shapes
:
Q
save/AssignVariableOp_4AssignVariableOpbeta_1save/Identity_4*
dtype0
Q
save/Identity_5Identitysave/RestoreV2:10*
T0*
_output_shapes
:
Q
save/AssignVariableOp_5AssignVariableOpbeta_2save/Identity_5*
dtype0
Q
save/Identity_6Identitysave/RestoreV2:11*
T0*
_output_shapes
:
P
save/AssignVariableOp_6AssignVariableOpdecaysave/Identity_6*
dtype0
Q
save/Identity_7Identitysave/RestoreV2:12*
T0*
_output_shapes
:
R
save/AssignVariableOp_7AssignVariableOpepsilonsave/Identity_7*
dtype0
`
save/Identity_8Identitysave/RestoreV2:13"/device:CPU:0*
_output_shapes
:*
T0	
^
save/AssignVariableOp_8AssignVariableOpitersave/Identity_8"/device:CPU:0*
dtype0	
Q
save/Identity_9Identitysave/RestoreV2:14*
_output_shapes
:*
T0
X
save/AssignVariableOp_9AssignVariableOplearning_ratesave/Identity_9*
dtype0
É
save/restore_all/NoOpNoOp^save/AssignVariableOp^save/AssignVariableOp_1^save/AssignVariableOp_2^save/AssignVariableOp_3^save/AssignVariableOp_4^save/AssignVariableOp_5^save/AssignVariableOp_6^save/AssignVariableOp_7^save/AssignVariableOp_9
^save/NoOp^save/NoOp_1^save/NoOp_2^save/NoOp_3^save/NoOp_4
H
save/restore_all/NoOp_1NoOp^save/AssignVariableOp_8"/device:CPU:0
J
save/restore_allNoOp^save/restore_all/NoOp^save/restore_all/NoOp_1
0
init_1NoOp^count_1/Assign^total_1/Assign"D
save/Const:0save/control_dependency:0save/restore_all 5 @F8"ň
trainable_variablesÚ×
x
dense/kernel:0dense/kernel/Assign"dense/kernel/Read/ReadVariableOp:0(2)dense/kernel/Initializer/random_uniform:08
g
dense/bias:0dense/bias/Assign dense/bias/Read/ReadVariableOp:0(2dense/bias/Initializer/zeros:08

dense_1/kernel:0dense_1/kernel/Assign$dense_1/kernel/Read/ReadVariableOp:0(2+dense_1/kernel/Initializer/random_uniform:08
o
dense_1/bias:0dense_1/bias/Assign"dense_1/bias/Read/ReadVariableOp:0(2 dense_1/bias/Initializer/zeros:08"Í
local_variablesšś
Y
	total_1:0total_1/Assigntotal_1/Read/ReadVariableOp:0(2total_1/Initializer/zeros:0
Y
	count_1:0count_1/Assigncount_1/Read/ReadVariableOp:0(2count_1/Initializer/zeros:0" 
cond_contextřô
ć

rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/cond_textrloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_t:0 *
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:0
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1
rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_t:0Ü
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Switch_1:1č
rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
Đr
tloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/cond_text_1rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_f:0*ú5
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Merge:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Merge:1
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch:1
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:1
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:0
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:2
§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1
¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim:0
 loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims:0
Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
Ťloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1
Śloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim:0
˘loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1:0
Ąloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat:0
Śloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims:0
Ľloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const:0
Ľloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Shape:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0
rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
sloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/switch_f:0
gloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank:0
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/rank:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch_1:0č
rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0rloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/pred_id:0
gloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/rank:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank/Switch:0
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:02,
,
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/cond_textloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0 *ç(
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:0
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:1
Žloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/DenseToDenseSetOperation:2
§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1
¤loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/dim:0
 loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims:0
Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
Ťloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1
Śloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/dim:0
˘loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1:0
Ąloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat/axis:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/concat:0
Śloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/num_invalid_dims:0
Ľloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Const:0
Ľloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like/Shape:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ones_like:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/x:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_t:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0Ťloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch_1:1Ř
Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims_1/Switch:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0Šloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch_1:1Ô
§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0§loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/has_invalid_dims/ExpandDims/Switch:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:02Ń
Î
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/cond_text_1loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0*

loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:1
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/switch_f:0
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/pred_id:0¤
loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/is_same_rank:0loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/has_valid_nonscalar_shape/Switch_1:0

oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/cond_textoloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t:0 *Ä
zloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency:0
oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0
ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_t:0â
oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0

qloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/cond_text_1oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f:0*ź
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch:0
wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1:0
wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2:0
wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_0:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_1:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_2:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_4:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_5:0
uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/data_7:0
|loss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/control_dependency_1:0
oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0
ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/switch_f:0
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0
ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge:0
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0â
oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0oloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/pred_id:0ä
iloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/weights/shape:0wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_1:0ă
hloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/values/shape:0wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_2:0ŕ
eloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_scalar:0wloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch_3:0é
ploss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/is_valid_shape/Merge:0uloss/dense_1_loss/mean_squared_error/weighted_loss/broadcast_weights/assert_broadcastable/AssertGuard/Assert/Switch:0"`
global_stepQO
M
iter:0iter/Assigniter/Read/ReadVariableOp:0(2iter/Initializer/zeros:0"Ž
	variables 
x
dense/kernel:0dense/kernel/Assign"dense/kernel/Read/ReadVariableOp:0(2)dense/kernel/Initializer/random_uniform:08
g
dense/bias:0dense/bias/Assign dense/bias/Read/ReadVariableOp:0(2dense/bias/Initializer/zeros:08

dense_1/kernel:0dense_1/kernel/Assign$dense_1/kernel/Read/ReadVariableOp:0(2+dense_1/kernel/Initializer/random_uniform:08
o
dense_1/bias:0dense_1/bias/Assign"dense_1/bias/Read/ReadVariableOp:0(2 dense_1/bias/Initializer/zeros:08
M
iter:0iter/Assigniter/Read/ReadVariableOp:0(2iter/Initializer/zeros:0
y
learning_rate:0learning_rate/Assign#learning_rate/Read/ReadVariableOp:0(2)learning_rate/Initializer/initial_value:0
Y
decay:0decay/Assigndecay/Read/ReadVariableOp:0(2!decay/Initializer/initial_value:0
]
beta_1:0beta_1/Assignbeta_1/Read/ReadVariableOp:0(2"beta_1/Initializer/initial_value:0
]
beta_2:0beta_2/Assignbeta_2/Read/ReadVariableOp:0(2"beta_2/Initializer/initial_value:0
a
	epsilon:0epsilon/Assignepsilon/Read/ReadVariableOp:0(2#epsilon/Initializer/initial_value:0*ă
evalÚ
3
dense_input$
dense_input:0˙˙˙˙˙˙˙˙˙

B
dense_1_target0
dense_1_target:0˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙˙9
metrics/accuracy/update_op
metric_op_wrapper:0 -
metrics/accuracy/value
Identity_11:0 ?
predictions/dense_1(
dense_1/BiasAdd:0˙˙˙˙˙˙˙˙˙
loss

loss/mul:0 tensorflow/supervised/eval*@
__saved_model_init_op'%
__saved_model_init_op
init_1
É
:
Add
x"T
y"T
z"T"
Ttype:
2	
B
AssignVariableOp
resource
value"dtype"
dtypetype
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
=
Mul
x"T
y"T
z"T"
Ttype:
2	

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
X
PlaceholderWithDefault
input"dtype
output"dtype"
dtypetype"
shapeshape
~
RandomUniform

shape"T
output"dtype"
seedint "
seed2int "
dtypetype:
2"
Ttype:
2	
@
ReadVariableOp
resource
value"dtype"
dtypetype
E
Relu
features"T
activations"T"
Ttype:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
:
Sub
x"T
y"T
z"T"
Ttype:
2	
q
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape
9
VarIsInitializedOp
resource
is_initialized
"serve*2.0.0-alpha02v1.12.0-9492-g2c319fb4158˛w
n
dense_inputPlaceholder*
dtype0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
*
shape:˙˙˙˙˙˙˙˙˙


-dense/kernel/Initializer/random_uniform/shapeConst*
valueB"
      *
_class
loc:@dense/kernel*
dtype0*
_output_shapes
:

+dense/kernel/Initializer/random_uniform/minConst*
valueB
 *č!ż*
_class
loc:@dense/kernel*
dtype0*
_output_shapes
: 

+dense/kernel/Initializer/random_uniform/maxConst*
valueB
 *č!?*
_class
loc:@dense/kernel*
dtype0*
_output_shapes
: 
Ě
5dense/kernel/Initializer/random_uniform/RandomUniformRandomUniform-dense/kernel/Initializer/random_uniform/shape*
T0*
_class
loc:@dense/kernel*
dtype0*
_output_shapes

:

Î
+dense/kernel/Initializer/random_uniform/subSub+dense/kernel/Initializer/random_uniform/max+dense/kernel/Initializer/random_uniform/min*
T0*
_class
loc:@dense/kernel*
_output_shapes
: 
ŕ
+dense/kernel/Initializer/random_uniform/mulMul5dense/kernel/Initializer/random_uniform/RandomUniform+dense/kernel/Initializer/random_uniform/sub*
_output_shapes

:
*
T0*
_class
loc:@dense/kernel
Ň
'dense/kernel/Initializer/random_uniformAdd+dense/kernel/Initializer/random_uniform/mul+dense/kernel/Initializer/random_uniform/min*
T0*
_class
loc:@dense/kernel*
_output_shapes

:


dense/kernelVarHandleOp*
shared_namedense/kernel*
_class
loc:@dense/kernel*
dtype0*
_output_shapes
: *
shape
:

i
-dense/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense/kernel*
_output_shapes
: 

dense/kernel/AssignAssignVariableOpdense/kernel'dense/kernel/Initializer/random_uniform*
_class
loc:@dense/kernel*
dtype0

 dense/kernel/Read/ReadVariableOpReadVariableOpdense/kernel*
_class
loc:@dense/kernel*
dtype0*
_output_shapes

:


dense/bias/Initializer/zerosConst*
valueB*    *
_class
loc:@dense/bias*
dtype0*
_output_shapes
:


dense/biasVarHandleOp*
dtype0*
_output_shapes
: *
shape:*
shared_name
dense/bias*
_class
loc:@dense/bias
e
+dense/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOp
dense/bias*
_output_shapes
: 
{
dense/bias/AssignAssignVariableOp
dense/biasdense/bias/Initializer/zeros*
_class
loc:@dense/bias*
dtype0

dense/bias/Read/ReadVariableOpReadVariableOp
dense/bias*
_class
loc:@dense/bias*
dtype0*
_output_shapes
:
h
dense/MatMul/ReadVariableOpReadVariableOpdense/kernel*
dtype0*
_output_shapes

:

r
dense/MatMulMatMuldense_inputdense/MatMul/ReadVariableOp*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
c
dense/BiasAdd/ReadVariableOpReadVariableOp
dense/bias*
dtype0*
_output_shapes
:
v
dense/BiasAddBiasAdddense/MatMuldense/BiasAdd/ReadVariableOp*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
S

dense/ReluReludense/BiasAdd*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
Ł
/dense_1/kernel/Initializer/random_uniform/shapeConst*
valueB"      *!
_class
loc:@dense_1/kernel*
dtype0*
_output_shapes
:

-dense_1/kernel/Initializer/random_uniform/minConst*
dtype0*
_output_shapes
: *
valueB
 *mż*!
_class
loc:@dense_1/kernel

-dense_1/kernel/Initializer/random_uniform/maxConst*
dtype0*
_output_shapes
: *
valueB
 *m?*!
_class
loc:@dense_1/kernel
Ň
7dense_1/kernel/Initializer/random_uniform/RandomUniformRandomUniform/dense_1/kernel/Initializer/random_uniform/shape*
dtype0*
_output_shapes

:*
T0*!
_class
loc:@dense_1/kernel
Ö
-dense_1/kernel/Initializer/random_uniform/subSub-dense_1/kernel/Initializer/random_uniform/max-dense_1/kernel/Initializer/random_uniform/min*
T0*!
_class
loc:@dense_1/kernel*
_output_shapes
: 
č
-dense_1/kernel/Initializer/random_uniform/mulMul7dense_1/kernel/Initializer/random_uniform/RandomUniform-dense_1/kernel/Initializer/random_uniform/sub*
T0*!
_class
loc:@dense_1/kernel*
_output_shapes

:
Ú
)dense_1/kernel/Initializer/random_uniformAdd-dense_1/kernel/Initializer/random_uniform/mul-dense_1/kernel/Initializer/random_uniform/min*
T0*!
_class
loc:@dense_1/kernel*
_output_shapes

:

dense_1/kernelVarHandleOp*
dtype0*
_output_shapes
: *
shape
:*
shared_namedense_1/kernel*!
_class
loc:@dense_1/kernel
m
/dense_1/kernel/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_1/kernel*
_output_shapes
: 

dense_1/kernel/AssignAssignVariableOpdense_1/kernel)dense_1/kernel/Initializer/random_uniform*!
_class
loc:@dense_1/kernel*
dtype0

"dense_1/kernel/Read/ReadVariableOpReadVariableOpdense_1/kernel*
dtype0*
_output_shapes

:*!
_class
loc:@dense_1/kernel

dense_1/bias/Initializer/zerosConst*
dtype0*
_output_shapes
:*
valueB*    *
_class
loc:@dense_1/bias

dense_1/biasVarHandleOp*
shape:*
shared_namedense_1/bias*
_class
loc:@dense_1/bias*
dtype0*
_output_shapes
: 
i
-dense_1/bias/IsInitialized/VarIsInitializedOpVarIsInitializedOpdense_1/bias*
_output_shapes
: 

dense_1/bias/AssignAssignVariableOpdense_1/biasdense_1/bias/Initializer/zeros*
_class
loc:@dense_1/bias*
dtype0

 dense_1/bias/Read/ReadVariableOpReadVariableOpdense_1/bias*
_class
loc:@dense_1/bias*
dtype0*
_output_shapes
:
l
dense_1/MatMul/ReadVariableOpReadVariableOpdense_1/kernel*
dtype0*
_output_shapes

:
u
dense_1/MatMulMatMul
dense/Reludense_1/MatMul/ReadVariableOp*
T0*'
_output_shapes
:˙˙˙˙˙˙˙˙˙
g
dense_1/BiasAdd/ReadVariableOpReadVariableOpdense_1/bias*
dtype0*
_output_shapes
:
|
dense_1/BiasAddBiasAdddense_1/MatMuldense_1/BiasAdd/ReadVariableOp*'
_output_shapes
:˙˙˙˙˙˙˙˙˙*
T0
,
predict/group_depsNoOp^dense_1/BiasAdd
U
ConstConst"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_1Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_2Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_3Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
\
Const_4Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB Bmodel
W
Const_5Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_6Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
W
Const_7Const"/device:CPU:0*
dtype0*
_output_shapes
: *
valueB B 
W
Const_8Const"/device:CPU:0*
valueB B *
dtype0*
_output_shapes
: 
¤
RestoreV2/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
r
RestoreV2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

	RestoreV2	RestoreV2Const_4RestoreV2/tensor_namesRestoreV2/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
B
IdentityIdentity	RestoreV2*
T0*
_output_shapes
:
I
AssignVariableOpAssignVariableOpdense/kernelIdentity*
dtype0
¤
RestoreV2_1/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_1/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_1	RestoreV2Const_4RestoreV2_1/tensor_namesRestoreV2_1/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_1IdentityRestoreV2_1*
T0*
_output_shapes
:
K
AssignVariableOp_1AssignVariableOp
dense/bias
Identity_1*
dtype0
Ś
RestoreV2_2/tensor_namesConst"/device:CPU:0*K
valueBB@B6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_2/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_2	RestoreV2Const_4RestoreV2_2/tensor_namesRestoreV2_2/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_2IdentityRestoreV2_2*
T0*
_output_shapes
:
O
AssignVariableOp_2AssignVariableOpdense_1/kernel
Identity_2*
dtype0
¤
RestoreV2_3/tensor_namesConst"/device:CPU:0*I
value@B>B4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
t
RestoreV2_3/shape_and_slicesConst"/device:CPU:0*
valueB
B *
dtype0*
_output_shapes
:

RestoreV2_3	RestoreV2Const_4RestoreV2_3/tensor_namesRestoreV2_3/shape_and_slices"/device:CPU:0*
_output_shapes
:*
dtypes
2
F

Identity_3IdentityRestoreV2_3*
T0*
_output_shapes
:
M
AssignVariableOp_3AssignVariableOpdense_1/bias
Identity_3*
dtype0
P
VarIsInitializedOpVarIsInitializedOpdense_1/kernel*
_output_shapes
: 
N
VarIsInitializedOp_1VarIsInitializedOp
dense/bias*
_output_shapes
: 
P
VarIsInitializedOp_2VarIsInitializedOpdense/kernel*
_output_shapes
: 
P
VarIsInitializedOp_3VarIsInitializedOpdense_1/bias*
_output_shapes
: 
d
initNoOp^dense/bias/Assign^dense/kernel/Assign^dense_1/bias/Assign^dense_1/kernel/Assign
Y
save/filename/inputConst*
valueB Bmodel*
dtype0*
_output_shapes
: 
n
save/filenamePlaceholderWithDefaultsave/filename/input*
shape: *
dtype0*
_output_shapes
: 
e

save/ConstPlaceholderWithDefaultsave/filename*
dtype0*
_output_shapes
: *
shape: 

save/Const_1Const*
dtype0*
_output_shapes
: *Í
valueĂBŔ Bš{"class_name": "Sequential", "config": {"layers": [{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "batch_input_shape": [null, 10], "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense", "trainable": true, "units": 5, "use_bias": true}}, {"class_name": "Dense", "config": {"activation": "linear", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense_1", "trainable": true, "units": 2, "use_bias": true}}], "name": "sequential"}}
Ö
save/Const_2Const*
valueB B{"class_name": "InputLayer", "config": {"batch_input_shape": [null, 10], "dtype": "float32", "name": "dense_input", "sparse": false}}*
dtype0*
_output_shapes
: 

save/Const_3Const*
dtype0*
_output_shapes
: *Ř
valueÎBË BÄ{"class_name": "Dense", "config": {"activation": "relu", "activity_regularizer": null, "batch_input_shape": [null, 10], "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense", "trainable": true, "units": 5, "use_bias": true}}
ř
save/Const_4Const*
dtype0*
_output_shapes
: *ť
valueąBŽ B§{"class_name": "Dense", "config": {"activation": "linear", "activity_regularizer": null, "bias_constraint": null, "bias_initializer": {"class_name": "Zeros", "config": {}}, "bias_regularizer": null, "dtype": "float32", "kernel_constraint": null, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "kernel_regularizer": null, "name": "dense_1", "trainable": true, "units": 2, "use_bias": true}}
ń
save/SaveV2/tensor_namesConst*¤
valueBB/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:
s
save/SaveV2/shape_and_slicesConst*
dtype0*
_output_shapes
:*#
valueBB B B B B B B B 
ą
save/SaveV2SaveV2
save/Constsave/SaveV2/tensor_namessave/SaveV2/shape_and_slicessave/Const_1save/Const_2save/Const_3dense/bias/Read/ReadVariableOp dense/kernel/Read/ReadVariableOpsave/Const_4 dense_1/bias/Read/ReadVariableOp"dense_1/kernel/Read/ReadVariableOp*
dtypes

2
}
save/control_dependencyIdentity
save/Const^save/SaveV2*
T0*
_class
loc:@save/Const*
_output_shapes
: 

save/RestoreV2/tensor_namesConst"/device:CPU:0*¤
valueBB/.ATTRIBUTES/OBJECT_CONFIG_JSONB&layer-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB3layer_with_weights-0/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB3layer_with_weights-1/.ATTRIBUTES/OBJECT_CONFIG_JSONB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
dtype0*
_output_shapes
:

save/RestoreV2/shape_and_slicesConst"/device:CPU:0*#
valueBB B B B B B B B *
dtype0*
_output_shapes
:
Â
save/RestoreV2	RestoreV2
save/Constsave/RestoreV2/tensor_namessave/RestoreV2/shape_and_slices"/device:CPU:0*
dtypes

2*4
_output_shapes"
 ::::::::

	save/NoOpNoOp

save/NoOp_1NoOp

save/NoOp_2NoOp
N
save/IdentityIdentitysave/RestoreV2:3*
_output_shapes
:*
T0
Q
save/AssignVariableOpAssignVariableOp
dense/biassave/Identity*
dtype0
P
save/Identity_1Identitysave/RestoreV2:4*
T0*
_output_shapes
:
W
save/AssignVariableOp_1AssignVariableOpdense/kernelsave/Identity_1*
dtype0

save/NoOp_3NoOp
P
save/Identity_2Identitysave/RestoreV2:6*
_output_shapes
:*
T0
W
save/AssignVariableOp_2AssignVariableOpdense_1/biassave/Identity_2*
dtype0
P
save/Identity_3Identitysave/RestoreV2:7*
T0*
_output_shapes
:
Y
save/AssignVariableOp_3AssignVariableOpdense_1/kernelsave/Identity_3*
dtype0
´
save/restore_allNoOp^save/AssignVariableOp^save/AssignVariableOp_1^save/AssignVariableOp_2^save/AssignVariableOp_3
^save/NoOp^save/NoOp_1^save/NoOp_2^save/NoOp_3

init_1NoOp"D
save/Const:0save/control_dependency:0save/restore_all 5 @F8"ň
trainable_variablesÚ×
x
dense/kernel:0dense/kernel/Assign"dense/kernel/Read/ReadVariableOp:0(2)dense/kernel/Initializer/random_uniform:08
g
dense/bias:0dense/bias/Assign dense/bias/Read/ReadVariableOp:0(2dense/bias/Initializer/zeros:08

dense_1/kernel:0dense_1/kernel/Assign$dense_1/kernel/Read/ReadVariableOp:0(2+dense_1/kernel/Initializer/random_uniform:08
o
dense_1/bias:0dense_1/bias/Assign"dense_1/bias/Read/ReadVariableOp:0(2 dense_1/bias/Initializer/zeros:08"č
	variablesÚ×
x
dense/kernel:0dense/kernel/Assign"dense/kernel/Read/ReadVariableOp:0(2)dense/kernel/Initializer/random_uniform:08
g
dense/bias:0dense/bias/Assign dense/bias/Read/ReadVariableOp:0(2dense/bias/Initializer/zeros:08

dense_1/kernel:0dense_1/kernel/Assign$dense_1/kernel/Read/ReadVariableOp:0(2+dense_1/kernel/Initializer/random_uniform:08
o
dense_1/bias:0dense_1/bias/Assign"dense_1/bias/Read/ReadVariableOp:0(2 dense_1/bias/Initializer/zeros:08*
serving_default
3
dense_input$
dense_input:0˙˙˙˙˙˙˙˙˙
3
dense_1(
dense_1/BiasAdd:0˙˙˙˙˙˙˙˙˙tensorflow/serving/predict*@
__saved_model_init_op'%
__saved_model_init_op
init_1