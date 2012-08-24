#ifndef __AppDelegate_H__
#define __AppDelegate_H__

#include "OgrePlatform.h"

#if OGRE_PLATFORM != OGRE_PLATFORM_APPLE
#error This header is for use with Mac OS X only
#endif

#ifdef __OBJC__

#import <Cocoa/Cocoa.h>

// All this does is suppress some messages in the run log.  NSApplication does not
// implement buttonPressed and apps without a NIB have no target for the action.
@implementation NSApplication (_suppressUnimplementedActionWarning)
- (void) buttonPressed:(id)sender
{
    /* Do nothing */
}
@end


#if defined(MAC_OS_X_VERSION_10_6) && MAC_OS_X_VERSION_MAX_ALLOWED >= MAC_OS_X_VERSION_10_6
@interface AppDelegate : NSObject <NSApplicationDelegate>
#else
@interface AppDelegate : NSObject
#endif
{
    NSTimer *mTimer;
    NSDate *mDate;
    double mLastFrameTime;
    double mStartTime;
}

- (void)startRendering;
- (void)renderOneFrame:(id)sender;

@property (retain) NSTimer *mTimer;
@property (nonatomic) double mLastFrameTime;
@property (nonatomic) double mStartTime;

@end

#if __LP64__
static id mAppDelegate;
#endif

@implementation AppDelegate

@synthesize mTimer;
@dynamic mLastFrameTime;
@dynamic mStartTime;

- (double)mLastFrameTime
{
    return mLastFrameTime;
}

- (void)setLastFrameTime:(double)frameInterval
{
    // Frame interval defines how many display frames must pass between each time the
    // display link fires. The display link will only fire 30 times a second when the
    // frame internal is two on a display that refreshes 60 times a second. The default
    // frame interval setting of one will fire 60 times a second when the display refreshes
    // at 60 times a second. A frame interval setting of less than one results in undefined
    // behavior.
    if (frameInterval >= 1)
    {
        mLastFrameTime = frameInterval;
    }
}

- (void)startRendering {
    
  NSAutoreleasePool * pool = [[NSAutoreleasePool alloc] init];
  mLastFrameTime = 1;
  mStartTime = 0;
  mTimer = nil;
   
  mTimer = [NSTimer scheduledTimerWithTimeInterval:(NSTimeInterval)(1.0f / 60.0f) * mLastFrameTime
  	    target:self
  	    selector:@selector(renderOneFrame:)
  	    userInfo:nil
  	    repeats:YES];
  //  mTimer = [NSTimer timerWithTimeInterval: (NSTimeInterval)(1.0f / 60.0f) * mLastFrameTime target: self selector: @selector(renderOneFrame:) userInfo: nil repeats: YES];
 
  //  [[NSRunLoop mainRunLoop] addTimer: mTimer forMode: NSRunLoopCommonModes];
  [pool release];
}
- (void)applicationDidFinishLaunching:(NSNotification *)application {
  mLastFrameTime = 1;
  mStartTime = 0;
  mTimer = nil;
    
  [self startRendering];
}


//- (NSApplicationTerminateReply)applicationShouldTerminate:(NSApplication*)sender 
//{
//  return NSTerminateCancel;
//}


- (void)renderOneFrame:(id)sender
{
  if(InputListener::getSingletonPtr()->viewerIsRunning() &&
     Ogre::Root::getSingletonPtr() && Ogre::Root::getSingleton().isInitialised())
    {
      mStartTime = InputListener::getSingletonPtr()->getTimer()->getMillisecondsCPU();
            
      InputListener::getSingletonPtr()->getKeyBoard()->capture();
      InputListener::getSingletonPtr()->getMouse()->capture();
            
      InputListener::getSingletonPtr()->updateViewer(mLastFrameTime);
      ViewerOgre3D::getSingleton().getOgreRoot()->renderOneFrame();
            
      mLastFrameTime = InputListener::getSingletonPtr()->getTimer()->getMillisecondsCPU() - mStartTime;
    }
  else
    {
         [mTimer invalidate];
         mTimer = nil;
         [NSApp performSelector:@selector(terminate:) withObject:nil afterDelay:0.0];
    }
}

- (void)dealloc {
    if(mTimer)
    {
        [mTimer invalidate];
        mTimer = nil;
    }
    
    [super dealloc];
}

@end

#endif

#endif
