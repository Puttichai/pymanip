## GLOBAL PARAMETERS

## DIRECTIONS (p means plus; m means minus)
pX = 0
pY = 1
pZ = 2
mX = 3
mY = 4
mZ = 5

GRASPDICT = {pX: "+X",
             pY: "+Y",
             pZ: "+Z",
             mX: "-X",
             mY: "-Y",
             mZ: "-Z"}


## preprocessing for rectangular object

def PossibleContactNormal(objectextents, dmax):
    objx = objectextents[0]
    objy = objectextents[1]
    objz = objectextents[2]
    normaldir = []
    
    if ((objx > dmax) and (objy > dmax) and (objz < dmax)):
        normaldir.append(pX)
        normaldir.append(pY)
        normaldir.append(mX)
        normaldir.append(mY)
    elif ((objx < dmax) and (objy > dmax) and (objz > dmax)):
        normaldir.append(pY)
        normaldir.append(pZ)
        normaldir.append(mY)
        normaldir.append(mZ)
    elif ((objx > dmax) and (objy < dmax) and (objz > dmax)):
        normaldir.append(pX)
        normaldir.append(pZ)
        normaldir.append(mX)
        normaldir.append(mZ)
    elif ((objx > dmax) and (objy < dmax) and (objz > dmax)):
        pass
    else:
        ## all six directions are possible
        normaldir.append(pX)
        normaldir.append(pY)
        normaldir.append(pZ)
        normaldir.append(mX)
        normaldir.append(mY)
        normaldir.append(mZ)

    return normaldir


def PossibleApproachingDirection(objectextents, dmax, gripperoffset):
    """
    PossibleApproachingDirection returns a set containing possible
    approaching direction of the gripper for each case of contact
    surfaces.
    
    return value: lib 

    lib[normaldir] is a set of possible approaching direction given
    that 'normaldir' is the normal direction of the surface in contact

    gripperoffset is to ensure that the object is high (w.r.t. the
    contact surface) enough such that the gripper can approach it from
    a side
    """
    objx = objectextents[0]
    objy = objectextents[1]
    objz = objectextents[2]

    lib = []

    ## contact normal is +X
    temp = []
    if ((objy < dmax) or (objz < dmax)):
        temp.append(pX)
    if (objx > gripperoffset):
        ## object's dimension allows the gripper to approach it from a
        ## side
        temp.append(pY)
        temp.append(mY)
        temp.append(pZ)
        temp.append(mZ)
    temp.sort()
    lib.append(temp)

    ## contact normal is +Y
    temp = []
    if ((objx < dmax) or (objz < dmax)):
        temp.append(pY)
    if (objy > gripperoffset):
        ## object's dimension allows the gripper to approach it from a
        ## side
        temp.append(pX)
        temp.append(mX)
        temp.append(pZ)
        temp.append(mZ)            
    temp.sort()
    lib.append(temp)
    
    ## contactnormal is +Z
    temp = []
    if ((objx < dmax) or (objy < dmax)):
        temp.append(pZ)
    if (objz > gripperoffset):
        ## object's dimension allows the gripper to approach it from a
        ## side
        temp.append(pX)
        temp.append(mX)
        temp.append(pY)
        temp.append(mY)
    temp.sort()
    lib.append(temp)

    ## contact normal is -X
    temp = []
    if ((objy < dmax) or (objz < dmax)):
        temp.append(mX)
    if (objx > gripperoffset):
        ## object's dimension allows the gripper to approach it from a
        ## side
        temp.append(pY)
        temp.append(mY)
        temp.append(pZ)
        temp.append(mZ)
    temp.sort()
    lib.append(temp)

    ## contact normal is -Y
    temp = []
    if ((objx < dmax) or (objz < dmax)):
        temp.append(mY)
    if (objy > gripperoffset):
        ## object's dimension allows the gripper to approach it from a
        ## side
        temp.append(pX)
        temp.append(mX)
        temp.append(pZ)
        temp.append(mZ)            
    temp.sort()
    lib.append(temp)
    
    ## contactnormal is -Z
    temp = []
    if ((objx < dmax) or (objy < dmax)):
        temp.append(mZ)
    if (objz > gripperoffset):
        ## object's dimension allows the gripper to approach it from a
        ## side
        temp.append(pX)
        temp.append(mX)
        temp.append(pY)
        temp.append(mY)
    temp.sort()
    lib.append(temp)

    return lib
    

def PossibleSlidingDirection(objectextents, dmax):
    """
    PossibleSlidingDirection returns a set containing possible sliding
    direction of the gripper for each case of approaching directions.
    
    return value: lib 

    lib[approachingdir] is a set of possible sliding direction given
    the 'approachingdir'.
    """
    objx = objectextents[0]
    objy = objectextents[1]
    objz = objectextents[2]

    lib = []
    
    ## approaching direction is +X
    temp = []
    if (objz < dmax):
        temp.append(pY)
        if (objy < dmax):
            temp.append(pZ)
    elif (objy < dmax):
        temp.append(pZ)
    lib.append(temp)

    ## approaching direction is +Y
    temp = []
    if (objz < dmax):
        temp.append(pX)
        if (objx < dmax):
            temp.append(pZ)
    elif (objx < dmax):
        temp.append(pZ)
    lib.append(temp)

    ## approaching direction is +Z
    temp = []
    if (objy < dmax):
        temp.append(pX)
        if (objx < dmax):
            temp.append(pY)
    elif (objx < dmax):
        temp.append(pY)
    lib.append(temp)

    ## approaching direction is -X
    lib.append(lib[pX])
    
    ## approaching direction is -Y
    lib.append(lib[pY])

    ## approaching direction is -Z
    lib.append(lib[pZ])
            
    return lib
            
