# Get selected animlayers
# Top one is the adjustment layer
# sum the curves from the lower layers

# if no selection, the top layer is the adjustment layer
# it will sum all layers below it

# It is assuming all layers involved with the selection are additive
# TODO: properly validate/calculate curves that are override

from maya import cmds, mel
import maya.api.OpenMaya as om
import maya.api.OpenMayaAnim as oma

DEBUG      = False
DO_SET     = True
CHANNELS   = [ 'translate'
             , 'rotate'
             , 'scale'
             ]
ATTRIBUTES = [ 'translateX'
             , 'translateY'
             , 'translateZ'
             , 'rotateX'
             , 'rotateY'
             , 'rotateZ'
             , 'scaleX'
             , 'scaleY'
             , 'scaleZ'
             ]
GRAPH_EDITOR = 'graphEditor1GraphEd'


if DEBUG:
    from pprint import pprint as pp

# Helper dict will create a new key if it doesn't already exist
class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value




def get_layers_to_process():
    ''' Returns a list of animation layers that will be processed
    '''

    all_layers = []
    root_layer = cmds.animLayer(query=True, root=True)
    if root_layer:
        all_layers.append(root_layer)
        all_layers.extend(cmds.animLayer(root_layer, q=True, children=True) or []) # top-down is last-to-first
    else:
        return None


    # Get selected animLayers
    selected_layers = cmds.treeView("AnimLayerTabanimLayerEditor", q=True, selectItem=True) or []

    layers_to_process = all_layers[:] # Copy the list
    if len(selected_layers) == 0 or selected_layers[-1] == 'BaseAnimation': # Why select base?? Treat it like selecting nothing, I guess.
        selected_layers = [layers_to_process.pop()] # Make it a list
    elif len(selected_layers) == 1:
        index = layers_to_process.index(selected_layers[-1])
        del layers_to_process[index:]
    elif len(selected_layers) > 1:
        layers_to_process = selected_layers[:-1]

    if DEBUG:
        print 'adjustment layer target is {0}. Summed layers are {1}'.format(selected_layers[-1], layers_to_process)

    adjustment_layer = selected_layers[-1]

    if not isinstance(layers_to_process, list): layers_to_process = [layers_to_process]
    # if not isinstance(adjustment_layer, list): adjustment_layer = [adjustment_layer]
    return adjustment_layer, layers_to_process


def get_attribute_curves(attribute, layer):
    ''' Returns dict
    {animationLayer : curveName}
    '''

    obj, attr = attribute.split('.')
    
    # Filter out non ATTRIBUTES
    if attr not in ATTRIBUTES:
        return None
    
    if attribute in cmds.animLayer(layer, q=True, attribute=True) or []:
        anim_curve = cmds.animLayer(layer, q=True, findCurveForPlug=attribute) or []
    
    if anim_curve:
        return {layer : anim_curve[0]}


def working_scale():
    # Maya ALWAYS works in CM. So scale internals by the working units in order to get the desired result!
    units = cmds.currentUnit(q=True, linear=True)
    if units == 'mm':
        return 0.1
    if units == 'cm':
        return 1
    if units == 'm':
        return 100
    if units == 'km':
        return 1000
    if units == 'in':
        return 2.54
    if units == 'ft':
        return 30.48
    if units == 'yd':
        return 91.44
    if units == 'mi':
        return 160934


def return_MFnAnimCurve(curve):
    msel = om.MSelectionList()
    msel.add(curve)
    mdep = msel.getDependNode(0)
    mcurve = oma.MFnAnimCurve(mdep)
    if not mcurve.name() == curve:
        return None
    return mcurve


def get_value_graph(mcurve, frange=None):
    if not frange:
        frange = get_curve_range(mcurve)

    values = []
    for frame in frange:
        # Get value at mtime, also feed it the current time uiUnit
        value = mcurve.evaluate(om.MTime(frame, om.MTime.uiUnit()))

        # Rotation curves are returned as Angular!
        if mcurve.animCurveType == oma.MFnAnimCurve.kAnimCurveTA:
            value = om.MAngle.internalToUI(value)
        else:
            # Translation curves are scaled by the WORKING UNITS! WTF?!
            value = value / working_scale()

        values.append(value)

    # Fighting precision errors :(
    rounded_values =  [round(x,10) for x in values]

    return rounded_values


def get_velocity_graph(values):
    velocity_graph = [0.0]
    for i in xrange(len(values)):
        if i > 0:
            current_value = values[i]
            previous_value = values[i-1]
            velocity_graph.append(abs(current_value - previous_value))
    return velocity_graph

def get_float_range(list_of_keys):

    float_range = [x * 1.0 for x in xrange(int(min(list_of_keys)), int(max(list_of_keys))+1)]

    set_range = set(float_range)
    set_range.update(list_of_keys)
    list_range = list(set_range)
    list_range.sort()
    
    return list_range
    
def get_curve_range(mcurve):
    key_times = []
    for key_index in xrange(int(mcurve.numKeys)):
        mtime = mcurve.input(key_index) # Get (time, timeType) at keyIndex
        key_times.append(mtime.value)

    int_range = [x * 1.0 for x in xrange(int(min(key_times)), int(max(key_times))+1)]

    set_range = set(int_range)
    set_range.update(key_times)
    list_range = list(set_range)
    list_range.sort()

    return list_range

def get_curve_ranges(mcurve):
    key_times = []
    for key_index in xrange(int(mcurve.numKeys)):
        mtime = mcurve.input(key_index) # Get (time, timeType) at keyIndex
        key_times.append(mtime.value)
    range_times = []
    for i, key in enumerate(key_times):
        if key != key_times[-1]:
            key_range = key_times[i], key_times[i+1]
            range_times.append(key_range)
    return range_times

def normalize_values(values, normal=100):
    # Returns a list of values where all values will add to 100
    # Normalize it to 100
    if abs(sum(values)) > 0.0:
        mult =  normal / abs(sum(values))
        return [(x * mult) for x in values]
    else:
        # cmds.error('Curve cannot be normalized. Values are flat.')
        return [0.0 for x in values]


# http://stackoverflow.com/q/3844948/
def is_equal(lst):
    return not lst or lst.count(lst[0]) == len(lst)


def map_from_to(x,a,b,c,d):
   y=(x-a)/(b-a)*(d-c)+c
   return y


def remap(old_value, old_min, old_max, new_min, new_max):
    old_range = (old_max - old_min)

    if old_range == 0:
        new_value = new_min
    else:
        new_range = (new_max - new_min)
        new_value = (((old_value - old_min) * new_range) / old_range) + new_min

    return new_value


def get_other_axis(attribute):
    ''' Takes a string,
        Replaces the X Y or Z with a list of the other two
    '''
    # Wow this is hacky wtf I'm sorry
    axis = ['X', 'Y', 'Z']
    attr = [x for x in axis if x in attribute]
    axis.remove(attr[0])
    return [attribute.replace(attr[0], axis[0]),
            attribute.replace(attr[0], axis[-1])]

def get_other_channel(channel):
    ''' Takes a string,
        Replaces the channel with a list of the other two
    '''
    # Wow this is hacky wtf I'm sorry
    channels = ['translate', 'rotate', 'scale']
    attr = [x for x in channels if x in channel]
    channels.remove(attr[0])
    return channels


def get_animated_attributes(node):
    import maya.OpenMaya as om1
    import maya.OpenMayaAnim as oma1

    # Get a MDagPath for the given node name:
    # node = 'pCube1'
    selList = om1.MSelectionList()
    selList.add(node)
    mDagPath = om1.MDagPath()
    selList.getDagPath(0, mDagPath)

    # Find all the animated attrs:
    mPlugArray = om1.MPlugArray()
    oma1.MAnimUtil.findAnimatedPlugs(mDagPath, mPlugArray)
    animCurves = []
    attribute_names = []

    # Find the curves ultimately connected to the attrs:
    for i in range(mPlugArray.length()):
        mPlugObj = om1.MPlug(mPlugArray[i])
        attribute_name = mPlugObj.name()
        attribute_names.append(attribute_name)

        # We could go on to capture all the animCurves of the plug
        # But this skips layer memberships.
        # We would need a better way traverse the connections.
        # Perhaps this holds the key...
        # https://discourse.techart.online/t/maya-animlayer-and-the-api/3510/4
        mObjArray = om1.MObjectArray()
        oma1.MAnimUtil.findAnimation(mPlugArray[i], mObjArray)
        for j in range(mObjArray.length()):
            depNodeFunc = om1.MFnDependencyNode(mObjArray[j])
            animCurves.append(depNodeFunc.name())

    # See what we found:
    # for ac in sorted(animCurves):
    #     print ac
    return sorted(attribute_names)

def apply_values(curve, values):
    # Do the magic, do the magic!
    for index, key in enumerate(values):
        cmds.keyframe(curve, index=(index,), valueChange=key, absolute=True)

def keywithmaxval(d): # https://stackoverflow.com/questions/268272/getting-key-with-maximum-value-in-dictionary
    """ a) create a list of the dict's keys and values;
        b) return the key with the max value"""
    v=list(d.values())
    k=list(d.keys())
    if is_equal(v):
        return None
    else:
        return k[v.index(max(v))]



    
def get_attribute_layer_curve(attribute, layer):
    """ Find the curve for the attribute on the specified layer
        If no curve is found, return the value of the attribute (assume unkeyed)
    """

    if not is_object_in_layer(attribute, layer):
        return None
    
    plug = None
    if layer == cmds.animLayer(q=True, root=True):
        # For the base animation layer, traverse the chain of animBlendNodes all
        # the way to the end.  The plug will be "inputA" on that last node.
        conns = cmds.listConnections(attribute, type='animBlendNodeBase', source=True, destination=False)
        blendNode = None
        while conns:
            blendNode = conns[0]
            conns = cmds.listConnections(blendNode, type='animBlendNodeBase', source=True, destination=False)
        plug = '{0}.inputA'.format(blendNode)
        return cmds.getAttr(plug)
    else:
        # For every layer other than the base animation layer, we can just use
        # the "animLayer" command.  Unfortunately the "layeredPlug" flag is
        # broken in Python in Maya 2016, so we have to use MEL.
        cmd = 'animLayer -q -layeredPlug "{0}" "{1}"'.format(attribute, layer)
        plug = mel.eval(cmd)
        return cmds.listConnections(plug)
    # return plug
        

def is_object_in_layer(obj, layer):
    """ Determine if the given object is in the given animation layer.
    """
    object_layer_members = cmds.animLayer([obj], q=True, affectedLayers=True) or []
    if layer in object_layer_members:
        return True
    return False
        




# ---------------------------------------------------------------------------- #
# Run commands

def run(smart=False, do_set=False):
    # We only process one layer and it's controls (or the selected controls in it)

    # Start by getting layers and which layer is the adjustment layer
    adjustment_layer, layers_to_process = get_layers_to_process() # Validates layer selection
    if not adjustment_layer:
        cmds.warning("No animation layers to process. Aborting!")
        # Eject if no layers
        return None

    adjustment_layer_members = cmds.animLayer( adjustment_layer
                                             , q=True
                                             , attribute=True)
    members = [x.split('.')[0] for x in adjustment_layer_members]
    members = list(set(members))
                                             
    if not adjustment_layer_members:
        cmds.warning("Adjustment layer {} has no members. Aborting!".format(adjustment_layer))
        # Eject if no layers
        return None

    objects = cmds.ls(sl=1)
    if not objects:
        # No selection? Fetch members of the adjustment layer
        cmds.warning("No object selected. Fetching members of {} instead.".format(adjustment_layer))
        objects = members
    else:
        if not bool(set(members) & set(objects)):
            print "No selected objects exist in the selected adjustment layer!"
            return None

    if not objects:
        # Still no objects? Abort.
        cmds.warning("{} contains no controls. Aborting!".format(adjustment_layer))
        return None

    # Eject non-member layers from layers_to_process
    layers_to_not_process = []
    for layer in layers_to_process:
        if layer == cmds.animLayer(q=True, root=True): continue
        layer_members = cmds.animLayer( layer
                                      , q=True
                                      , attribute=True)
        members = [x.split('.')[0] for x in layer_members]
        members = list(set(members))
        if not bool(set(members) & set(objects)):
            print "Didn't find {0} in {1}.".format(members, layer)
            layers_to_not_process.append(layer)
    
    for layer in layers_to_not_process:
        layers_to_process.remove(layer)

                    

    # Validation done

    # ======================================================================= #
    # BEGIN
    # ======================================================================= #
    
    adjustment_keys = set()
    ctrl_curves_to_process = Vividict()
    
    # adjustment_layer = 'AnimLayer4'
    # adjustment_layer = 'BaseAnimation'

    # This section will populate the dictionary like so:
    # control 
    # ∟ attribute
    #   ∟ animation layer
    #     ∟ animation curve OR static value

    for obj in objects:
        
        animated_attributes = get_animated_attributes(obj) # The API version
        
        for layer in layers_to_process + [adjustment_layer]:
            
            for attribute in animated_attributes: 
                
                obj, attr = attribute.split('.')

                if attr not in ATTRIBUTES:
                    continue # Whitelisting attributes for now
                
                # if 'rotate' not in attr: continue
                
                if attribute in adjustment_layer_members:
                    if layer == cmds.animLayer(q=True, root=True): # BaseAnimation is treated differently... thanks Maya
                        # Now we traverse the tree going from the top animLayer down to the base
                        connections = cmds.listConnections(attribute, type='animBlendNodeBase', source=True, destination=False)
                        blend_node = None
                        while connections:
                            blend_node = connections[0]
                            connections = cmds.listConnections(blend_node, type='animBlendNodeBase', source=True, destination=False)
                        plug = '{0}.inputA'.format(blend_node) # We hit base
                        # print 'Plug is {}'.format(plug)
                        curve = cmds.listConnections(plug) or []
                        if curve:
                            ctrl_curves_to_process[obj][attr][layer] = curve
                            # print 'Curve is {}'.format(curve)
                        else:
                            # Treat rotations differently on the base as well... thanks maya
                            if cmds.nodeType(blend_node) == 'animBlendNodeAdditiveRotation':
                                if 'X' in attr:
                                    ctrl_curves_to_process[obj][attr][layer] = cmds.getAttr(plug)[0][0]
                                    # print cmds.getAttr(plug)[0][0]
                                if 'Y' in attr:
                                    ctrl_curves_to_process[obj][attr][layer] = cmds.getAttr(plug)[0][1]
                                    # print cmds.getAttr(plug)[0][1]
                                if 'Z' in attr:
                                    ctrl_curves_to_process[obj][attr][layer] = cmds.getAttr(plug)[0][2]
                                    # print cmds.getAttr(plug)[0][2]
                                
                            else: # Everything else is fine
                                ctrl_curves_to_process[obj][attr][layer] = cmds.getAttr(plug)
                    else:
                        plug = cmds.animLayer(layer, q=True, layeredPlug=attribute)
                        # print 'Plug is {}'.format(plug)
                        curve = cmds.animLayer(layer, q=True, findCurveForPlug=attribute)
                        if curve:
                            # print 'Curve is {}'.format(curve)
                            ctrl_curves_to_process[obj][attr][layer] = curve
                        else:
                            # print cmds.getAttr(plug.replace('.inputB', '.inputA'))
                            try:
                                ctrl_curves_to_process[obj][attr][layer] = cmds.getAttr(plug.replace('.inputB', '.inputA'))
                            except: 
                                print obj, attr, layer

                        if layer == adjustment_layer:
                        
                            if curve:
                                ctrl_curves_to_process[obj][attr][layer] = curve
                                # Now get in there and get all the keys on the adjustment layer curves
                                # Do this per-object so we can preserve 'the pose' on those frames
                                keyframes = cmds.keyframe(curve, q=True)
                                for key in keyframes:
                                    adjustment_keys.add(key)
                            else:
                                cmds.error("No adjustment curve found for {}".format(attribute))
                                continue # Something is fucked if this happens

        if not ctrl_curves_to_process[obj]: continue # Eject ghosts

    if not ctrl_curves_to_process: return False # How does this happen?

    # At this point, we have the curve names of objects on the 
    # adjustment layer, and keys of all the objects on this layer.
    # So we can composite adjustment ranges between these keys
    
    adjustment_key_ranges = []
    adjustment_keys_sorted = sorted(adjustment_keys)
    for index, key in enumerate(adjustment_keys_sorted):
        if not index == len(adjustment_keys_sorted) - 1:
            adjustment_key_ranges.append([key, adjustment_keys_sorted[index+1]])
    
    # Working calculation range
    if not adjustment_keys:
        cmds.error("Could not find any adjustment keys on {}".format(adjustment_layer))
        return False

    calculation_range = get_float_range(adjustment_keys)
    
    # Now we need to calculate the layers_to_process between these ranges
    # We can query all curves between these ranges to get value graphs
    # If we don't find a curve (no key on BaseAnimation for example), we can grab the flat value
    
    for obj in ctrl_curves_to_process.keys():
        
        # get_rotates    = False
        # get_translates = False
        
        for attr in ctrl_curves_to_process[obj].keys():
            
            # if smart:
            #     if 'rotate' in attr:
            #         get_rotates = True
            #     if 'translate' in attr:
            #         get_translates = True

            for layer, destination in ctrl_curves_to_process[obj][attr].items():

                if not isinstance(destination, list): 
                    destination = ctrl_curves_to_process[obj][attr][layer] = [destination]

                if isinstance(destination[0], float):
                    float_range = [destination[0] for x in calculation_range]

                elif isinstance(destination[0], unicode):
                    float_range = []
                    for time in calculation_range:
                        value = cmds.keyframe(destination[0], q=True, valueChange=True, eval=True, time=(time,))[0]
                        float_range.append(value)

                else:
                    cmds.error("Something went horribly wrong with {0}, {1}.".format(layer, destination))
                    continue

                if layer == adjustment_layer and is_equal(float_range):
                    if DEBUG:
                        print "Ejecting flat adjustment curve for {0}.{1}".format(obj, attr)
                    del ctrl_curves_to_process[obj][attr][layer]
                    continue
                    
                ctrl_curves_to_process[obj][attr][layer].append(float_range)
    
    
    # ======================================================================= #
    # Begin calculation of the curve data

    value_graphs = Vividict()
    
    for obj in ctrl_curves_to_process.keys():

        for attr in ctrl_curves_to_process[obj].keys():
            
            composite_velocity_graphs = []
            
            for layer, destination in ctrl_curves_to_process[obj][attr].items():
                
                if isinstance(destination[0], float):
                    continue

                elif isinstance(destination[0], unicode):
                    api_curve = return_MFnAnimCurve(destination[0])
                    value_graph = get_value_graph(api_curve, calculation_range)

                    if is_equal(value_graph):
                        continue
  
                    if layer == adjustment_layer:
                        value_graphs[obj][attr]['adjustment_graph'] = value_graph
                        value_graphs[obj][attr]['adjustment_curve'] = destination[0]
                    else:
                        composite_velocity_graphs.append(get_velocity_graph(value_graph))
                    # print value_graph
            
            if composite_velocity_graphs:
                for graph in composite_velocity_graphs:
                    for i, value in enumerate(composite_velocity_graphs):
                        if i != 0:
                            for x,_ in enumerate(value):
                                composite_velocity_graphs[0][x] += value[x]
                value_graphs[obj][attr]['composite_graph'] = composite_velocity_graphs[0]

    #  len(calculation_range)
    adjustment_range = range(int(adjustment_key_ranges[0][0]), int(adjustment_key_ranges[-1][-1])+1)
    
    for obj in value_graphs.keys():
        for attr in value_graphs[obj].keys():

            adjustment_curve = value_graphs[obj][attr]['adjustment_curve']
            adjustment_graph = value_graphs[obj][attr]['adjustment_graph']
            composite_graph  = value_graphs[obj][attr]['composite_graph']
            
            if not composite_graph or not adjustment_curve or not adjustment_graph:
                # DO SMART SHIT
                if DEBUG:
                    print "Ejecting {}. No values to calculate.".format(attr)
                continue

            new_value_curve = []
            frame_march = []
            
            for frange in adjustment_key_ranges:

                frame_range = range(int(frange[0]), int(frange[1])+1)

                normalized_velocity_graph = normalize_values(composite_graph[calculation_range.index(frange[0]):calculation_range.index(frange[1])+1])
                
                # if is_equal(normalized_velocity_graph): continue # How did this end up here?

                sum_percentage = 0.0

                for index, value in enumerate(frame_range):
                    sum_percentage += normalized_velocity_graph[index]
                    new_value = map_from_to(sum_percentage, 0, 100, adjustment_graph[calculation_range.index(frange[0])], adjustment_graph[calculation_range.index(frange[1])])
                    if value not in frame_march:
                        new_value_curve.append(new_value)
                        frame_march.append(value) # I do this to skip the repeat frames between sets - those keys already exist anyway

            # Now set the keys
            # Do the magic, DO THE MAGIC!
            if do_set:
                if DEBUG:
                    print "Running adjustment on {}.".format(adjustment_curve)
                for index, time in enumerate(adjustment_range):
                    cmds.setKeyframe(adjustment_curve, animLayer=adjustment_layer, time=(time,), value=new_value_curve[index])


    '''
    for obj in ctrl_curves_to_process.keys():
        total_adjustment_range = []
        for attribute, layers in ctrl_curves_to_process[obj].items():
            adjustment_curve  = layers[adjustment_layer][0]
            adjustment_range  = get_curve_range(return_MFnAnimCurve(adjustment_curve))
            total_adjustment_range.extend(adjustment_range)
        ctrl_curves_to_process[obj]['total_adjustment_range'] = sorted(list(set(total_adjustment_range)))


    for obj in ctrl_curves_to_process.keys():
        attributes_to_delete = []
        total_adjustment_range = ctrl_curves_to_process[obj]['total_adjustment_range']
        del ctrl_curves_to_process[obj]['total_adjustment_range'] # this key was just a passenger. No longer needed.
        for attribute, layers in ctrl_curves_to_process[obj].items():

            layer_composite = []

            adjustment_curve  = layers[adjustment_layer][0]
            adjustment_range  = get_curve_range(return_MFnAnimCurve(adjustment_curve))
            adjustment_ranges = get_curve_ranges(return_MFnAnimCurve(adjustment_curve))

            for layer, curve in layers.items():

                api_curve = return_MFnAnimCurve(curve[0])

                value_graph = get_value_graph(api_curve, total_adjustment_range)

                if layer == adjustment_layer:
                    if is_equal(value_graph):
                        # If the adjustment layer is flat, nothing can be done here.
                        attributes_to_delete.append(attribute)
                        # continue

                    ctrl_curves_to_process[obj][attribute]['adjustment_curve']  = curve[0]
                    ctrl_curves_to_process[obj][attribute]['adjustment_layer']  = layer
                    ctrl_curves_to_process[obj][attribute]['adjustment_range']  = adjustment_range
                    ctrl_curves_to_process[obj][attribute]['adjustment_ranges'] = adjustment_ranges
                    ctrl_curves_to_process[obj][attribute]['adjustment_values'] = value_graph

                else:
                    velocity_graph = get_velocity_graph(value_graph)
                    layer_composite.append(velocity_graph)

            # Composite the layers together
            for i, value in enumerate(layer_composite):
                if i != 0:
                    for x,_ in enumerate(value):
                        layer_composite[0][x] += value[x]

            ctrl_curves_to_process[obj][attribute]['composite_velocity'] = layer_composite[0]

            # Normalize it for final consumption
            # normalized_velocity_graph = normalize_values(layer_composite[0])
            # ctrl_curves_to_process[obj][attribute]['composite_velocity_normalized'] = normalized_velocity_graph

        # Skip the attributes that cannot be 'adjusted'
        for attr in attributes_to_delete:
            del ctrl_curves_to_process[obj][attr]


    print "Running operation on \n{}".format('\n'.join(ctrl_curves_to_process[obj].keys()))
    for obj in ctrl_curves_to_process.keys():
        for attribute, data in ctrl_curves_to_process[obj].items():
            adjustment_curve = data['adjustment_curve']
            adjustment_layer = data['adjustment_layer']
            adjustment_range = data['adjustment_range']
            adjustment_ranges = data['adjustment_ranges']
            adjustment_values = data['adjustment_values']
            composite_velocity_graph  = data['composite_velocity']

            if not adjustment_ranges: continue # Nothing to adjust

            if is_equal(composite_velocity_graph):
                if not smart:
                    print "Aint smart with the {}".format(attribute)
                    continue
                # We do the clever shit here
                # Get neighboring axis and composite them
                redundant_keys = {}
                sum_keys = {}
                for attr in ATTRIBUTES:
                    if attr in attribute:
                        current_axis = attr
                        other_axis = get_other_axis(current_axis)
                        for axis in other_axis:
                            new_axis = attribute.replace(current_axis, axis)
                            composite_velocity_graph = ctrl_curves_to_process[obj][new_axis]['composite_velocity']

                            if sum(composite_velocity_graph) == 0: continue # worthless flat curve
                            sum_keys[new_axis] = sum(composite_velocity_graph)
                            redundants = 0.0
                            for index, value in enumerate(composite_velocity_graph):
                                if index == 0: continue
                                if value == composite_velocity_graph[index - 1]:
                                    redundants += 1
                            redundant_keys[new_axis] = redundants

                if redundant_keys: # Here we check to see if the velocity ever goes flat. This is important in animation
                    key = keywithmaxval(redundant_keys)
                    if not key:
                        key = keywithmaxval(sum_keys)
                    if not key:
                        key = redundant_keys.keys()[0]
                    composite_velocity_graph = ctrl_curves_to_process[obj][key]['composite_velocity']

            if is_equal(composite_velocity_graph):
                redundant_keys = {}
                for channel in CHANNELS:
                    if channel in attribute:
                        current_channel = channel
                        other_channels = get_other_channel(channel)
                        for other_channel in other_channels:
                            new_channel = attribute.replace(current_channel, other_channel)
                            for attr in ATTRIBUTES:
                                if attr in new_channel:
                                    current_axis = attr
                                    other_axis = get_other_axis(current_axis)
                                    other_axis.append(current_axis)
                                    for axis in other_axis:
                                        new_axis = new_channel.replace(current_axis, axis)
                                        composite_velocity_graph = ctrl_curves_to_process[obj][new_axis]['composite_velocity']

                                        if sum(composite_velocity_graph) == 0:
                                            continue # worthless flat curve

                                        redundants = 0.0
                                        for index, value in enumerate(composite_velocity_graph):
                                            if index == 0: continue
                                            if value == composite_velocity_graph[index - 1]:
                                                redundants += 1
                                        redundant_keys[new_axis] = redundants

                if redundant_keys:
                    key = keywithmaxval(redundant_keys)
                    if not key:
                        key = redundant_keys.keys()[0]
                    composite_velocity_graph = ctrl_curves_to_process[obj][key]['composite_velocity']

            if is_equal(composite_velocity_graph):
                continue # Nothing can be done in this attribute group

                # continue # SKIP IT for now - we don't have the clever shit installed

            new_value_curve = []
            frame_march = []
            for frange in adjustment_ranges:

                frame_range = range(int(frange[0]), int(frange[1])+1)

                normalized_velocity_graph = normalize_values(composite_velocity_graph[adjustment_range.index(frange[0]):adjustment_range.index(frange[1])+1])

                sum_percentage = 0.0

                for index, value in enumerate(frame_range):
                    sum_percentage += normalized_velocity_graph[index]
                    new_value = map_from_to(sum_percentage, 0, 100, adjustment_values[adjustment_range.index(frange[0])], adjustment_values[adjustment_range.index(frange[1])])
                    if value not in frame_march:
                        new_value_curve.append(new_value)
                        frame_march.append(value) # I do this to skip the repeat frames between sets - those keys already exist anyway

            # Now set the keys
            # Do the magic, DO THE MAGIC!
            if DO_SET:
                for index, time in enumerate(adjustment_range):
                    cmds.setKeyframe(adjustment_curve, animLayer=adjustment_layer, time=(time,), value=new_value_curve[index])
    '''

# ---------------------------------------------------------------------------- #
# Bunch of dev shit here

def num_reversals(values):
    reverals = []
    begin = False
    falling = False

    for index, value in enumerate(values):
        if index == 0: # ignore first key
            continue

        if value == values[index-1]: # ignore redunant keys
            continue

        # First direction change
        if begin == False:
            if value < values[index-1]:
                reverals.append(values[index-1])
                falling = True
            elif value > values[index-1]:
                reverals.append(values[index-1])
                falling = False
            begin = True
            continue

        if value < values[index-1] and falling == False:
            reverals.append(values[index-1])
            falling = True
            # continue
        elif value > values[index-1] and falling == True:
            reverals.append(values[index-1])
            falling = False
        continue
    return reverals


def get_peaks_valleys(curve, frange=None):
    if isinstance(curve, str):
        mcurve = return_MFnAnimCurve(curve)
    elif isinstance(curve, oma.MFnAnimCurve):
        mcurve = curve
    else:
        cmds.error("Could not fetch curve from {}".format(curve))
        return None

    if not frange:
        frange = get_curve_range(mcurve)

    frame_difference = frange[-1] - frange[0]
    frame_difference = 1 if frame_difference == 0 else frame_difference

    value_graph = get_value_graph(mcurve)
    value_graph_times = []
    for index in range(mcurve.numKeys):
        time = mcurve.input(index)
        value_graph_times.append(time.value)

    # Skewing to right to match left value
    value_graph_skewed = skew_curve(curve)


def skew_values(values):
    frame_difference = len(values) - 1
    frame_difference = 1 if frame_difference == 0 else frame_difference

    offset_value = values[-1] - values[0] # The difference from first to last frame

    value_graph_skewed = []
    for index, value in enumerate(values):
        # frame = frange[index]

        time_slope = 1 - ((index - 1) / frame_difference) # Count from 1.0 to 0.0
        pivot_value = value - offset_value
        # Basically, just multiply it by the offset then multiply THAT by how far down the frange we are
        new_value = ((value - pivot_value) * time_slope) + pivot_value

        value_graph_skewed.append(new_value)

    return value_graph_skewed


def skew_curve(curve, frange=None):
    if isinstance(curve, str):
        mcurve = return_MFnAnimCurve(curve)
    elif isinstance(curve, oma.MFnAnimCurve):
        mcurve = curve
    else:
        cmds.error("Could not fetch curve from {}".format(curve))
        return None

    if not frange:
        frange = get_curve_range(mcurve)

    frame_difference = frange[-1] - frange[0]
    frame_difference = 1 if frame_difference == 0 else frame_difference

    value_graph = get_value_graph(mcurve)

    first_value = mcurve.value(0)
    last_value = mcurve.value(mcurve.numKeys - 1)
    # offset_value = first_value - last_value # The difference from first to last frame
    offset_value = last_value - first_value # The difference from first to last frame

    value_graph_skewed = []
    for index, value in enumerate(value_graph):
        frame = frange[index]

        time_slope = 1 - ((frame - frange[0]) / frame_difference) # Count from 1.0 to 0.0
        pivot_value = value - offset_value
        # Basically, just multiply it by the offset then multiply THAT by how far down the frange we are
        new_value = ((value - pivot_value) * time_slope) + pivot_value

        value_graph_skewed.append(new_value)

    return value_graph_skewed

def get_curve_intensity(curve):
    if isinstance(curve, str):
        mcurve = return_MFnAnimCurve(curve)
    elif isinstance(curve, oma.MFnAnimCurve):
        mcurve = curve
    else:
        cmds.error("Could not fetch curve from {}".format(curve))
        return None

    curve_data = {}

    value_graph = get_value_graph(mcurve)
    velocity_graph = get_velocity_graph(value_graph)

    value_graph_skewed = skew_curve(curve)
    reversals = num_reversals(value_graph_skewed)

    pivot_value = value_graph_skewed[0]
    peaks = []
    valleys = []
    for point in reversals:
        if point > pivot_value:
            peaks.append(point)
        elif point < pivot_value:
            valleys.append(point)

    redundants = 0.0
    for index, value in enumerate(velocity_graph):
        if index == 0: continue
        if value == velocity_graph[index - 1]:
            redundants += 1


    # draw a straight line from beginning to end
    # Every time you get a reversal on the top side, it is a peak
    num_peaks = len(peaks)
    num_valleys = len(valleys)
    # how big are the peaks vs valleys?
    highest_value = max(peaks)
    lowest_value = min(valleys)

    # hottest moment?
    highest_velocity = max(velocity_graph)
    total_change = sum(velocity_graph)

    # roll it into a data set
    curve_data['redundants']       = redundants
    # curve_data['num_peaks']      = num_peaks
    # curve_data['total_change']     = total_change
    # curve_data['num_valleys']    = num_valleys
    # curve_data['lowest_value']   = lowest_value
    # curve_data['num_reversals']  = len(reversals)
    # curve_data['highest_value']  = highest_value
    curve_data['highest_velocity'] = highest_velocity


    return curve_data


def compare_curve_intensities(curve1, curve2):
    # Counts the number of signals data1 beats over data2
    # Returns the winning curve
    data1 = get_curve_intensity(curve1)
    data2 = get_curve_intensity(curve2)

    data1_winnings = []
    for k in data1.keys():
        data1_winner = data1[k] > data2[k]
        data1_winnings.append(data1_winner)
    if data1_winnings.count(True) > data1_winnings.count(False):
        return curve1
    else:
        return curve2

# foo = compare_curve_intensities(curve1 = 'pCube1_rotateZ', curve2 = 'pCube1_rotateX')

def get_selected_curves():
    # get the key selection
    if not cmds.animCurveEditor(GRAPH_EDITOR, exists=True):
        cmds.error("No GraphEditor found.")
        return # Cannot find graph editor?

    if not cmds.animCurveEditor(GRAPH_EDITOR, q=True, areCurvesSelected=True):
        cmds.warning("Must select some keys to fit.")
        return

    selected_curves = cmds.keyframe(q=True, selected=True, name=True) or []

    return selected_curves


def get_curve_data():
    curves = get_selected_curves()
    anim_data = {}
    all_frames = []
    for curve in curves:
        selected_frames = cmds.keyframe(curve, q=True, selected=True, timeChange=True)
        all_frames.extend(selected_frames)

        # selected_index = cmds.keyframe(curve, q=True, selected=True, indexValue=True)
        selected_values = cmds.keyframe(curve, q=True, selected=True, valueChange=True)
        anim_data[curve] = [selected_frames, selected_values]

    first_frame = min(all_frames)
    last_frame = max(all_frames)


# ---------------------------------------------------------------------------- #
# Developer section

if __name__ == '__main__':
    run(smart=False, do_set=True)
    # pass
