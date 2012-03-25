
#import <Cocoa/Cocoa.h>

extern int nbgl_main_apple(int argc, const char* argv[]);

int main(int argc, const char* argv[])
{
#if 0
    [[NSAutoreleasePool alloc] init];
    [NSApplication sharedApplication];

    /*
    [[[NSNib alloc] initWithContentsOfURL:[NSURL URLWithString:@"MainMenu.nib"]] instantiateNibWithOwner:NSApp topLevelObjects:nil];
    */

    ProcessSerialNumber psn = { 0, kCurrentProcess };
    TransformProcessType(&psn, kProcessTransformToForegroundApplication);
    SetFrontProcess(&psn);

    [NSApp finishLaunching];
#endif

    return nbgl_main_apple(argc, argv);
}

